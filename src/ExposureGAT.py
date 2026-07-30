
"""
ExposureGAT.py

Training loop for ExposureGAT.

Objective
---------
The normal -> tumor transition (delta_h, delta_alpha) should be larger
for patients with higher smoking exposure (pack_years).

This is a pairwise ranking objective:
  - For every pair (i, j) in a batch where pack_years[i] > pack_years[j] + tau,
    enforce that transition_score[i] > transition_score[j].

Transition score per patient:
  node_contribution : mean L2-norm of delta_h   per graph
  edge_contribution : mean |delta_alpha|          per graph
  combined          : node_weight * node_score + edge_weight * edge_score

Loss:
  mean softplus(trans[j] - trans[i] + margin)
  over all valid pairs (i, j) in the batch.

Model selection:
  Best model = highest Spearman correlation between transition score
  and pack_years on the validation set.

Outputs saved to out_dir/
  best_model.pth      state dict of best model
  last_model.pth      state dict after final epoch
  training_log.csv    per-epoch metrics
  run_config.json     all hyperparameters
  
%run ExposureGAT \
  --data_dir GNN_Input_Data \
  --out_dir  GNN_Output\
  --epochs   100 
  """


import os
import copy
import json
import random
import argparse
import traceback
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Subset
from torch_geometric.data import Batch

from BuildGraph import MultiOmicsDataset
from PrepareDataset import ExposureGAT


# =============================================================================
# Reproducibility
# =============================================================================
def set_seed(seed: int = 42) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


# =============================================================================
# Collate function
# =============================================================================
def collate_fn(batch):
    normals, tumors, exposures, pids = zip(*batch)
    batch_normal = Batch.from_data_list(list(normals))
    batch_tumor  = Batch.from_data_list(list(tumors))
    exposure     = torch.stack(list(exposures), dim=0)   # [B, exp_dim]
    return batch_normal, batch_tumor, exposure, list(pids)


# =============================================================================
# Exposure maps for controls
# =============================================================================
class ExposureOverrideDataset:
    """
    Thin wrapper around MultiOmicsDataset that replaces the exposure tensor
    returned by __getitem__ with a fixed patient-level exposure mapping.

    Used for shuffled-exposure control. Molecular graphs stay unchanged; only
    the patient -> supplied exposure assignment is changed.
    """

    def __init__(self, base_dataset: MultiOmicsDataset, exposure_map: dict):
        self.base_dataset = base_dataset
        self.exposure_map = {
            str(k): np.asarray(v, dtype=np.float32) for k, v in exposure_map.items()
        }

    def __len__(self) -> int:
        return len(self.base_dataset)

    def __getitem__(self, idx: int):
        normal_graph, tumor_graph, exposure, pid = self.base_dataset[idx]
        pid = str(pid)
        if pid not in self.exposure_map:
            raise KeyError(f"No exposure override found for patient {pid}")

        raw_exp = self.exposure_map[pid]
        exposure = torch.tensor(
            (raw_exp - self.base_dataset._exp_mean) / self.base_dataset._exp_std,
            dtype=torch.float32,
        )
        return normal_graph, tumor_graph, exposure, pid


def _patient_raw_exposure_table(dataset: MultiOmicsDataset) -> pd.DataFrame:
    pid_col = dataset.patient_id_col
    exp_cols = dataset.exposure_cols
    patient_meta = dataset.meta.drop_duplicates(subset=[pid_col]).copy()
    patient_meta[pid_col] = patient_meta[pid_col].astype(str)
    return patient_meta.set_index(pid_col)[exp_cols].astype(np.float32)


def make_true_standardized_exposure_map(dataset: MultiOmicsDataset) -> Dict[str, np.ndarray]:
    """Patient -> TRUE standardized exposure using scalers fitted on train set."""
    raw_table = _patient_raw_exposure_table(dataset)
    true_map: Dict[str, np.ndarray] = {}
    for pid, row in raw_table.iterrows():
        raw = row.to_numpy(dtype=np.float32)
        true_map[str(pid)] = ((raw - dataset._exp_mean) / dataset._exp_std).astype(np.float32)
    return true_map


def make_shuffled_exposure_map(
    dataset: MultiOmicsDataset,
    train_pids: List[str],
    val_pids: List[str],
    seed: int,
) -> Tuple[dict, pd.DataFrame]:
    """
    Create a fixed patient -> shuffled RAW exposure mapping for a negative control.

    Shuffling is done separately within train and validation splits. This keeps
    the exposure distribution in each split but destroys patient-exposure pairing.
    """
    rng = np.random.RandomState(seed)
    exp_cols = dataset.exposure_cols
    raw_table = _patient_raw_exposure_table(dataset)

    exposure_map = {}
    records = []

    def _shuffle_group(pids: List[str], split_name: str) -> None:
        pids = [str(p) for p in pids]
        raw = raw_table.loc[pids, exp_cols].to_numpy(dtype=np.float32)
        shuffled = raw.copy()
        rng.shuffle(shuffled)  # shuffle rows, preserving multicolumn exposure vectors

        for pid, true_exp, shuf_exp in zip(pids, raw, shuffled):
            exposure_map[pid] = shuf_exp
            rec = {"patient_id": pid, "split": split_name}
            for j, col in enumerate(exp_cols):
                rec[f"true_{col}"] = float(true_exp[j])
                rec[f"shuffled_{col}"] = float(shuf_exp[j])
            records.append(rec)

    _shuffle_group(train_pids, "train")
    _shuffle_group(val_pids, "val")
    return exposure_map, pd.DataFrame(records)


def true_exposure_tensor(
    pids: List[str],
    true_exp_map: Dict[str, np.ndarray],
    device: torch.device,
) -> torch.Tensor:
    vals = [true_exp_map[str(pid)] for pid in pids]
    return torch.tensor(np.stack(vals, axis=0), dtype=torch.float32, device=device)


# =============================================================================
# Transition score
# =============================================================================
def compute_transition_score(
    delta_h:     torch.Tensor,
    delta_alpha: torch.Tensor,
    edge_index:  torch.Tensor,
    node_batch:  torch.Tensor,
    n_graphs:    int,
    node_weight: float = 1.0,
    edge_weight: float = 1.0,
) -> torch.Tensor:
    """Compute one scalar normal->tumor transition score per patient graph."""
    node_mag   = delta_h.norm(dim=1)
    node_score = torch.zeros(n_graphs, device=delta_h.device)
    node_cnt   = torch.zeros(n_graphs, device=delta_h.device)
    node_score.scatter_add_(0, node_batch, node_mag)
    node_cnt.scatter_add_(0, node_batch, torch.ones_like(node_mag))
    node_score = node_score / node_cnt.clamp_min(1.0)

    src        = edge_index[0]
    edge_graph = node_batch[src]
    edge_mag   = delta_alpha.abs()
    edge_score = torch.zeros(n_graphs, device=delta_h.device)
    edge_cnt   = torch.zeros(n_graphs, device=delta_h.device)
    edge_score.scatter_add_(0, edge_graph, edge_mag)
    edge_cnt.scatter_add_(0, edge_graph, torch.ones_like(edge_mag))
    edge_score = edge_score / edge_cnt.clamp_min(1.0)

    return node_weight * node_score + edge_weight * edge_score


# =============================================================================
# Pairwise ranking loss
# =============================================================================
def pairwise_ranking_loss(
    trans:    torch.Tensor,
    exposure: torch.Tensor,
    tau:      float = 1.0,
    margin:   float = 0.01,
) -> Tuple[torch.Tensor, int]:
    """
    For every ordered pair (i, j) where exposure[i] > exposure[j] + tau,
    penalise if trans[i] <= trans[j].

    Vectorised — O(n^2) in memory but avoids O(n^2) Python loop overhead.
    For batch_size=16 the matrix is 16x16=256 elements, negligible.
    Safe up to batch_size ~512 before memory becomes a concern.

    Previous version used a Python double for-loop which was correct
    but slow and incompatible with torch.compile contexts.
    """
    smk = exposure[:, 0]                                    # [B]

    # [B,1] - [1,B] = [B,B]  entry (i,j) = smk[i] - smk[j]
    smk_diff = smk.unsqueeze(1) - smk.unsqueeze(0)          # [B,B]

    # mask[i,j] = True  ->  trans[i] should rank above trans[j]
    mask = smk_diff > tau                                    # [B,B] bool

    n_pairs = int(mask.sum().item())
    if n_pairs == 0:
        return torch.tensor(0.0, device=trans.device), 0

    # entry (i,j) = trans[j] - trans[i]  (positive when i fails to beat j)
    trans_diff = trans.unsqueeze(0) - trans.unsqueeze(1)    # [B,B]

    loss = F.softplus(trans_diff + margin)[mask].mean()
    return loss, n_pairs


def _safe_spearman(x: List[float], y: List[float]) -> Tuple[float, float]:
    if len(x) < 4 or len(y) < 4:
        return float("nan"), float("nan")
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan"), float("nan")
    r, p = spearmanr(x, y)
    return float(r), float(p)


# =============================================================================
# One epoch (train or val)
# =============================================================================
def run_epoch(
    model:        torch.nn.Module,
    loader:       DataLoader,
    optimizer:    torch.optim.Optimizer,
    args:         argparse.Namespace,
    device:       torch.device,
    train:        bool,
    true_exp_map: Dict[str, np.ndarray],
) -> dict:

    model.train(train)
    ctx = torch.enable_grad() if train else torch.no_grad()

    total_loss       = 0.0
    all_trans        = []
    all_target_smk   = []  # exposure used for ranking loss/model selection
    all_true_smk     = []  # original true smoking, even in shuffled control
    all_supplied_smk = []  # exposure returned by dataloader; shuffled if control
    all_da_mean      = []
    n_batches        = 0
    n_skipped        = 0

    with ctx:
        for batch_normal, batch_tumor, supplied_exposure, pids in loader:
            try:
                batch_normal      = batch_normal.to(device)
                batch_tumor       = batch_tumor.to(device)
                supplied_exposure = supplied_exposure.to(device)
                true_exposure     = true_exposure_tensor(pids, true_exp_map, device)

                # Choose what the model sees and what the ranking loss uses.
                if args.no_exposure_model:
                    model_exposure  = torch.zeros_like(true_exposure)
                    target_exposure = true_exposure
                else:
                    model_exposure  = supplied_exposure
                    target_exposure = supplied_exposure

                for name, t in [
                    ("normal.x", batch_normal.x),
                    ("tumor.x", batch_tumor.x),
                    ("model_exposure", model_exposure),
                    ("target_exposure", target_exposure),
                    ("true_exposure", true_exposure),
                ]:
                    if not torch.isfinite(t).all():
                        raise ValueError(f"Non-finite values in {name}")

                out = model(batch_normal, batch_tumor, model_exposure)

                trans = compute_transition_score(
                    delta_h     = out["delta_h"],
                    delta_alpha = out["delta_alpha"],
                    edge_index  = out["edge_index"],
                    node_batch  = batch_tumor.batch,
                    n_graphs    = batch_tumor.num_graphs,
                    node_weight = args.node_weight,
                    edge_weight = args.edge_weight,
                )

                loss, n_pairs = pairwise_ranking_loss(
                    trans    = trans,
                    exposure = target_exposure,
                    tau      = args.pair_tau,
                    margin   = args.pair_margin,
                )

                # Always collect metrics, even if no valid ranking pairs.
                all_trans.extend(trans.detach().cpu().numpy().tolist())
                all_target_smk.extend(target_exposure[:, 0].detach().cpu().numpy().tolist())
                all_true_smk.extend(true_exposure[:, 0].detach().cpu().numpy().tolist())
                all_supplied_smk.extend(supplied_exposure[:, 0].detach().cpu().numpy().tolist())
                all_da_mean.append(float(out["delta_alpha"].abs().mean().detach().cpu()))

                if n_pairs == 0:
                    n_skipped += 1
                    if train:
                        continue
                    n_batches += 1
                    continue

                if not torch.isfinite(loss):
                    raise ValueError("Non-finite loss")

                if train:
                    optimizer.zero_grad(set_to_none=True)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
                    optimizer.step()

                total_loss += float(loss.detach().cpu())
                n_batches += 1

            except Exception as e:
                n_skipped += 1
                print(f"  [{'train' if train else 'val'}] skipped batch: {e}")
                traceback.print_exc()
                continue

    n = max(n_batches, 1)

    target_r, target_p     = _safe_spearman(all_trans, all_target_smk)
    true_r, true_p         = _safe_spearman(all_trans, all_true_smk)
    supplied_r, supplied_p = _safe_spearman(all_trans, all_supplied_smk)

    return {
        "loss":               total_loss / n,
        # Backward-compatible names: model selection uses target exposure.
        "spearman_r":         target_r,
        "spearman_p":         target_p,
        "spearman_target_r":  target_r,
        "spearman_target_p":  target_p,
        "spearman_true_r":    true_r,
        "spearman_true_p":    true_p,
        "spearman_supplied_r": supplied_r,
        "spearman_supplied_p": supplied_p,
        "mean_trans":         float(np.nanmean(all_trans))   if all_trans   else float("nan"),
        "std_trans":          float(np.nanstd(all_trans))    if all_trans   else float("nan"),
        "mean_da":            float(np.nanmean(all_da_mean)) if all_da_mean else float("nan"),
        "n_batches":          n_batches,
        "n_skipped":          n_skipped,
    }


# =============================================================================
# Stratified train / val split by true smoking exposure
# =============================================================================
def stratified_patient_split(
    dataset:  MultiOmicsDataset,
    val_frac: float,
    seed:     int,
) -> Tuple[List[str], List[str]]:
    exposure_col = dataset.exposure_cols[0]
    meta = dataset.meta.drop_duplicates(subset=[dataset.patient_id_col]).copy()
    pids = meta[dataset.patient_id_col].astype(str).values
    smk  = meta[exposure_col].fillna(0).values

    try:
        bins = pd.qcut(smk, q=3, labels=False, duplicates="drop")
    except ValueError:
        bins = np.zeros(len(smk), dtype=int)

    rng = np.random.RandomState(seed)
    train_pids, val_pids = [], []

    for b in np.unique(bins):
        stratum = pids[bins == b].copy()
        rng.shuffle(stratum)
        n_val = max(1, int(np.ceil(val_frac * len(stratum))))
        val_pids.extend(stratum[:n_val].tolist())
        train_pids.extend(stratum[n_val:].tolist())

    return train_pids, val_pids


# =============================================================================
# Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="ExposureGAT training with controls")

    # Paths
    ap.add_argument("--data_dir", type=str, default="GNN_Input_Data")
    ap.add_argument("--out_dir",  type=str, default="GNN_Output")
    ap.add_argument("--seed",     type=int, default=42)

    # Dataset scaling floors
    ap.add_argument("--meth_min_scale", type=float, default=0.005)
    ap.add_argument("--cnv_min_scale",  type=float, default=0.05)

    # Training
    ap.add_argument("--epochs",       type=int,   default=100)
    ap.add_argument("--batch_size",   type=int,   default=16)
    ap.add_argument("--lr",           type=float, default=1e-3)
    ap.add_argument("--weight_decay", type=float, default=1e-4)
    ap.add_argument("--val_frac",     type=float, default=0.20)
    ap.add_argument("--grad_clip",    type=float, default=5.0)

    # Ranking loss
    ap.add_argument("--pair_tau",    type=float, default=1.0)
    ap.add_argument("--pair_margin", type=float, default=0.01)
    ap.add_argument("--node_weight", type=float, default=1.0)
    ap.add_argument("--edge_weight", type=float, default=1.0)

    # LR scheduler
    ap.add_argument("--sched_patience", type=int,   default=8)
    ap.add_argument("--sched_factor",   type=float, default=0.5)
    ap.add_argument("--min_lr",         type=float, default=1e-6)

    # Model architecture
    ap.add_argument("--hidden_dim",  type=int,   default=64)
    ap.add_argument("--out_dim",     type=int,   default=32)
    ap.add_argument("--heads",       type=int,   default=4)
    ap.add_argument("--dropout",     type=float, default=0.1)
    ap.add_argument("--tanh_scale",  type=float, default=0.5)
    ap.add_argument("--bias_hidden", type=int,   default=32)
    ap.add_argument("--nudge_scale", type=float, default=0.1,
                    help="Degree nudge strength. 0.0 disables degree nudge.")

    # Controls
    ap.add_argument("--shuffle_exposure", action="store_true",
                    help="Control 1: train with shuffled supplied exposure. "
                         "Validation reports Spearman vs both shuffled target and TRUE exposure.")
    ap.add_argument("--shuffle_exposure_seed", type=int, default=42)
    ap.add_argument("--no_exposure_model", action="store_true",
                    help="Control 3: feed zero exposure to model, but train/evaluate ranking against TRUE exposure.")

    args = ap.parse_args()

    if args.shuffle_exposure and args.no_exposure_model:
        raise ValueError("Use either --shuffle_exposure or --no_exposure_model, not both.")

    set_seed(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # ------------------------------------------------------------------
    # Dataset and split
    # ------------------------------------------------------------------
    dataset = MultiOmicsDataset(
        data_dir       = args.data_dir,
        exposure_cols  = ("smoking_intensity",),
        meth_min_scale = args.meth_min_scale,
        cnv_min_scale  = args.cnv_min_scale,
        verbose        = True,
    )

    train_pids, val_pids = stratified_patient_split(dataset, args.val_frac, args.seed)
    print(f"Patients — total: {len(dataset.patient_ids)}  "
          f"train: {len(train_pids)}  val: {len(val_pids)}")

    smk_vals = dataset.meta.drop_duplicates(
        subset=[dataset.patient_id_col]
    )[dataset.exposure_cols[0]].fillna(0)
    print(f"Exposure ({dataset.exposure_cols[0]}) — "
          f"min={smk_vals.min():.1f}  median={smk_vals.median():.1f}  "
          f"max={smk_vals.max():.1f}  "
          f"patients with >0: {(smk_vals > 0).sum()}/{len(smk_vals)}")

    dataset.fit_scalers(train_patient_ids=train_pids)
    true_exp_map = make_true_standardized_exposure_map(dataset)

    train_idx = [i for i, p in enumerate(dataset.patient_ids) if p in set(train_pids)]
    val_idx   = [i for i, p in enumerate(dataset.patient_ids) if p in set(val_pids)]

    dataset_for_loader = dataset
    if args.shuffle_exposure:
        exposure_map, shuffle_df = make_shuffled_exposure_map(
            dataset=dataset,
            train_pids=train_pids,
            val_pids=val_pids,
            seed=args.shuffle_exposure_seed,
        )
        shuffle_path = os.path.join(args.out_dir, "shuffled_exposure_map.csv")
        shuffle_df.to_csv(shuffle_path, index=False)
        dataset_for_loader = ExposureOverrideDataset(dataset, exposure_map)
        print("*** CONTROL 1: shuffled supplied exposure ***")
        print("Model input and ranking target use shuffled exposure.")
        print("Metrics also report Spearman vs TRUE exposure.")
        print(f"Shuffled exposure map saved to: {shuffle_path}")
    elif args.no_exposure_model:
        print("*** CONTROL 3: no-exposure graph-only model ***")
        print("Model receives zero exposure; ranking target is TRUE exposure.")
    else:
        print("*** REAL EXPOSURE-GATED MODEL ***")

    train_loader = DataLoader(
        Subset(dataset_for_loader, train_idx),
        batch_size = args.batch_size,
        shuffle    = True,
        collate_fn = collate_fn,
        drop_last  = True,
    )

    val_loader = DataLoader(
        Subset(dataset_for_loader, val_idx),
        batch_size = max(1, len(val_idx)),
        shuffle    = False,
        collate_fn = collate_fn,
    )

    # ------------------------------------------------------------------
    # Model
    # ------------------------------------------------------------------
    model = ExposureGAT(
        cpg_in      = 1,
        gene_in     = 3,
        exp_dim     = len(dataset.exposure_cols),
        hidden_dim  = args.hidden_dim,
        out_dim     = args.out_dim,
        heads       = args.heads,
        dropout     = args.dropout,
        tanh_scale  = args.tanh_scale,
        bias_hidden = args.bias_hidden,
        nudge_scale = args.nudge_scale,
        n_nodes              = dataset.n_nodes,
        edge_index_for_degree= dataset.edge_index,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model parameters: {n_params:,}")

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr           = args.lr,
        weight_decay = args.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode     = "max",
        factor   = args.sched_factor,
        patience = args.sched_patience,
        min_lr   = args.min_lr,
    )

    # ------------------------------------------------------------------
    # Training loop
    # ------------------------------------------------------------------
    best_spearman = -float("inf")
    best_epoch    = -1
    logs: List[dict] = []

    best_model_path = os.path.join(args.out_dir, "best_model.pth")
    last_model_path = os.path.join(args.out_dir, "last_model.pth")
    log_path        = os.path.join(args.out_dir, "training_log.csv")
    cfg_path        = os.path.join(args.out_dir, "run_config.json")

    for epoch in range(1, args.epochs + 1):
        tr = run_epoch(model, train_loader, optimizer, args, device, train=True,
                       true_exp_map=true_exp_map)
        va = run_epoch(model, val_loader, optimizer, args, device, train=False,
                       true_exp_map=true_exp_map)

        current_lr = optimizer.param_groups[0]["lr"]
        scheduler.step(va["spearman_target_r"] if np.isfinite(va["spearman_target_r"]) else -1.0)

        row = {"epoch": epoch, "lr": current_lr}
        for prefix, metrics in [("train", tr), ("val", va)]:
            for k, v in metrics.items():
                row[f"{prefix}_{k}"] = v
        logs.append(row)
        pd.DataFrame(logs).to_csv(log_path, index=False)

        vt = va["spearman_target_r"]
        vv = va["spearman_true_r"]
        trt = tr["spearman_target_r"]
        print(
            f"Epoch {epoch:03d} | "
            f"TrainLoss {tr['loss']:.4f} | ValLoss {va['loss']:.4f} | "
            f"Val Spearman target {vt:.4f} (p={va['spearman_target_p']:.3f}) | "
            f"Val Spearman TRUE {vv:.4f} (p={va['spearman_true_p']:.3f}) | "
            f"Train Spearman target {trt:.4f} | "
            f"Val trans μ={va['mean_trans']:.3f} σ={va['std_trans']:.3f} | "
            f"Val |Δα| {va['mean_da']:.4f} | "
            f"Batches tr={tr['n_batches']} skip={tr['n_skipped']} "
            f"val={va['n_batches']} skip={va['n_skipped']} | "
            f"LR {current_lr:.2e}"
        )

        if np.isfinite(vt) and vt > best_spearman:
            best_spearman = vt
            best_epoch    = epoch
            torch.save(copy.deepcopy(model.state_dict()), best_model_path)
            print(f"  -> Best model saved (epoch {epoch}, target SpearmanR={best_spearman:.4f})")
            # Issue 2 fix: warn clearly when the saved model is a control.
            # The shuffled-exposure best model ranks patients by shuffled labels,
            # NOT by true smoking biology.  Loading it in Step_05 for biological
            # analysis would produce meaningless results.
            if args.shuffle_exposure:
                print("     *** CONTROL MODEL — ranked by SHUFFLED exposure. "
                      "Do NOT use best_model.pth for biological analysis. ***")
            elif args.no_exposure_model:
                print("     *** CONTROL MODEL — no exposure input. "
                      "Transition scores reflect graph structure only. ***")

    torch.save(model.state_dict(), last_model_path)

    # ------------------------------------------------------------------
    # Issue 4 fix: persist architecture constants that are NOT argparse args
    # but ARE required by Step_05 to reconstruct the model correctly.
    #
    # vars(args) only captures command-line flags.  n_nodes is derived from
    # dataset topology and nudge_scale is already in args, but n_nodes is
    # not — if Step_05 uses default n_nodes=0 it will raise a ValueError.
    # ------------------------------------------------------------------
    config_dict = vars(args).copy()
    config_dict["n_nodes"]  = dataset.n_nodes
    config_dict["n_params"] = n_params
    with open(cfg_path, "w") as f:
        json.dump(config_dict, f, indent=2)

    print("\n" + "=" * 60)
    print("Training complete.")
    print(f"Best epoch    : {best_epoch}  (target SpearmanR={best_spearman:.4f})")
    print(f"Best model    : {best_model_path}")
    print(f"Log           : {log_path}")
    print(f"Config        : {cfg_path}")
    print("=" * 60)

    if hasattr(model, "report_degree_thresholds"):
        print("\nLearned degree thresholds:")
        model.report_degree_thresholds()

    # ------------------------------------------------------------------
    # End-of-run summary banner — control status unmissable in logs
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    if args.shuffle_exposure:
        best_row = pd.DataFrame(logs).iloc[pd.DataFrame(logs)["val_spearman_target_r"].idxmax()]
        print(f"Best val Spearman TARGET: {best_row['val_spearman_target_r']:.4f}")
        print(f"Best val Spearman TRUE:   {best_row['val_spearman_true_r']:.4f}")
        print("CONTROL RUN — SHUFFLED EXPOSURE")
        print("best_model.pth: ranked by SHUFFLED labels.")
        print("Use spearman_true_r in training_log.csv to assess")
        print("whether graph structure alone carries smoking signal.")
        print("DO NOT use best_model.pth for biological analysis.")
    elif args.no_exposure_model:
        print("CONTROL RUN — NO EXPOSURE MODEL")
        print("best_model.pth: trained WITHOUT smoking exposure input.")
        print("Transition scores reflect molecular graph changes only.")
        print("DO NOT use best_model.pth as the main analysis model.")
    else:
        print("REAL MODEL RUN — exposure-gated, full architecture.")
        print("best_model.pth: safe for Step_05 biological analysis.")
    print("=" * 60)


if __name__ == "__main__":
    main()

