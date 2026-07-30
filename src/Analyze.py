#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze.py

Clean post-training result extraction for ExposureGAT.


It does:
  1. Run inference for each patient.
  2. Compute patient-level transition scores.
  3. Compute node transition magnitude: NodeDelta = mean ||delta_h||.
  4. Compute edge rewiring magnitude: EdgeDelta = mean |delta_alpha|.
  5. Build a high-rewiring network from top EdgeDelta edges.
  6. Compute node network support in that high-rewiring network:
       - RewiredDegree
       - WeightedRewiredDegree
       - LocalClusteringCoefficient
  7. Compute one final priority node score:
       FinalNodePriorityScore / NetworkWeightedNodeScore =
           w_delta   * scaled(NodeDelta)
         + w_degree  * scaled(WeightedRewiredDegree)
         + w_smoking * scaled(NodeSmokingAssociationScore)

     IMPORTANT: this final score uses scaled magnitudes, not ranks.
     NodeSmokingAssociationScore is computed as absolute Spearman correlation
     between patient-specific node transition magnitude and the selected
     smoking/exposure column.
  8. Export simple reporting/pathway-ready tables.
"""

import os
import json
import argparse
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

import torch
from torch_geometric.data import Batch

from BuildGraph import MultiOmicsDataset
from PrepareDataset import ExposureGAT


# =============================================================================
# Helpers
# =============================================================================


def rank_norm(x: np.ndarray) -> np.ndarray:
    """Rank-normalize a numeric vector to [0, 1]. Kept only for diagnostics."""
    x = np.asarray(x, dtype=float)
    out = np.zeros(len(x), dtype=float)
    mask = np.isfinite(x)
    if mask.sum() == 0:
        return out
    r = rankdata(x[mask], method="average")
    if len(r) == 1:
        out[mask] = 1.0
    else:
        out[mask] = (r - 1.0) / (len(r) - 1.0)
    return out


def minmax01(x: np.ndarray) -> np.ndarray:
    """Min-max scale a numeric vector to [0, 1], preserving score magnitude."""
    x = np.asarray(x, dtype=float)
    out = np.zeros(len(x), dtype=float)
    mask = np.isfinite(x)
    if mask.sum() == 0:
        return out

    xmin = np.nanmin(x[mask])
    xmax = np.nanmax(x[mask])
    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-12:
        return out

    out[mask] = (x[mask] - xmin) / (xmax - xmin)
    return out


def compute_node_smoking_association(
    node_mag: np.ndarray,
    exposures: np.ndarray,
    exposure_cols: List[str],
    preferred_exposure_col: str | None = None,
) -> Tuple[np.ndarray, str]:
    """
    Compute absolute Spearman association between patient-specific node transition
    magnitude and the selected smoking/exposure column.

    node_mag shape: patients x nodes.
    exposures shape: patients x exposure columns.
    """
    if exposures is None or np.asarray(exposures).size == 0 or len(exposure_cols) == 0:
        return np.zeros(node_mag.shape[1], dtype=float), "none"

    exposure_cols_str = [str(c) for c in exposure_cols]

    if preferred_exposure_col is not None and preferred_exposure_col in exposure_cols_str:
        exposure_idx = exposure_cols_str.index(preferred_exposure_col)
        selected_col = preferred_exposure_col
    else:
        preferred_names = [
            "smoking_intensity",
            "smoking",
            "pack_years",
            "pack_years_smoked",
            "cigarettes_per_day",
            "exposure",
            "Exposure",
        ]
        exposure_idx = 0
        selected_col = exposure_cols_str[0]
        for name in preferred_names:
            if name in exposure_cols_str:
                exposure_idx = exposure_cols_str.index(name)
                selected_col = name
                break

    y = np.asarray(exposures[:, exposure_idx], dtype=float)
    scores = np.zeros(node_mag.shape[1], dtype=float)

    for j in range(node_mag.shape[1]):
        x = np.asarray(node_mag[:, j], dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 4 or np.nanstd(x[ok]) < 1e-12 or np.nanstd(y[ok]) < 1e-12:
            scores[j] = 0.0
            continue
        r, _ = spearmanr(x[ok], y[ok])
        scores[j] = abs(float(r)) if np.isfinite(r) else 0.0

    return scores, selected_col


def clean_gene_label(label: str) -> str:
    """Convert GENE:1234 to 1234 and CpG:cg... to cg... for simple exports."""
    label = str(label)
    if label.startswith("GENE:"):
        return label.split(":", 1)[1]
    if label.startswith("CpG:") or label.startswith("CPG:"):
        return label.split(":", 1)[1]
    return label


def node_type_from_index(node_idx: int, n_genes: int) -> str:
    return "gene" if int(node_idx) < int(n_genes) else "CpG"


def edge_type_name(edge_type_value: int) -> str:
    return "PPI" if int(edge_type_value) == 0 else "CpG-Gene"


def local_clustering_coefficients(n_nodes: int, undirected_pairs: List[Tuple[int, int]]) -> np.ndarray:
    """
    Compute local clustering coefficients for an undirected graph.
    Nodes with degree < 2 receive 0.
    Implemented without networkx/igraph to keep Step 05 lightweight.
    """
    adj = [set() for _ in range(n_nodes)]
    edge_set = set()

    for u, v in undirected_pairs:
        u = int(u)
        v = int(v)
        if u == v:
            continue
        a, b = sorted((u, v))
        edge_set.add((a, b))
        adj[a].add(b)
        adj[b].add(a)

    clustering = np.zeros(n_nodes, dtype=float)

    for n in range(n_nodes):
        neighbors = list(adj[n])
        k = len(neighbors)
        if k < 2:
            clustering[n] = 0.0
            continue

        possible = k * (k - 1) / 2.0
        observed = 0
        for i in range(k):
            for j in range(i + 1, k):
                a, b = sorted((neighbors[i], neighbors[j]))
                if (a, b) in edge_set:
                    observed += 1
        clustering[n] = observed / possible if possible > 0 else 0.0

    return clustering


def write_csv(df: pd.DataFrame, out_dir: str, filename: str) -> None:
    path = os.path.join(out_dir, filename)
    df.to_csv(path, index=False)
    print(f"  Saved {filename} ({len(df)} rows)")


# =============================================================================
# Inference
# =============================================================================


def run_inference(model: ExposureGAT, dataset: LungMultiOmicsDataset, device: torch.device) -> Dict:
    """Run per-patient inference and collect delta_h, delta_alpha, alpha_n, alpha_t."""
    model.eval()
    model.to(device)

    single_normal = dataset[0][0]
    edge_index_np = single_normal.edge_index.cpu().numpy()
    edge_type_np = single_normal.edge_type.cpu().numpy()

    all_pids: List[str] = []
    all_exposures: List[np.ndarray] = []
    all_delta_h: List[np.ndarray] = []
    all_delta_alpha: List[np.ndarray] = []
    all_alpha_n: List[np.ndarray] = []
    all_alpha_t: List[np.ndarray] = []

    print(f"  Running per-patient inference (n={len(dataset)} patients)...")

    with torch.no_grad():
        for idx in range(len(dataset)):
            normal_g, tumor_g, exposure, pid = dataset[idx]
            normal_b = Batch.from_data_list([normal_g]).to(device)
            tumor_b = Batch.from_data_list([tumor_g]).to(device)
            exp_t = exposure.unsqueeze(0).to(device)

            out = model(normal_b, tumor_b, exp_t)

            all_pids.append(pid)
            all_exposures.append(exposure.cpu().numpy())
            all_delta_h.append(out["delta_h"].cpu().numpy())
            all_delta_alpha.append(out["delta_alpha"].cpu().numpy())
            all_alpha_n.append(out["alpha_n"].cpu().numpy())
            all_alpha_t.append(out["alpha_t"].cpu().numpy())

            if (idx + 1) % 20 == 0:
                print(f"    {idx + 1}/{len(dataset)} patients done")

    return {
        "patient_ids": all_pids,
        "exposures": np.stack(all_exposures, axis=0),
        "delta_h": all_delta_h,
        "delta_alpha": all_delta_alpha,
        "alpha_n": all_alpha_n,
        "alpha_t": all_alpha_t,
        "edge_index": edge_index_np,
        "edge_type": edge_type_np,
        "n_nodes": dataset.n_nodes,
        "n_genes": dataset.n_genes,
    }


# =============================================================================
# Core tables
# =============================================================================


def make_patient_table(raw: Dict, exposure_cols: List[str], node_weight: float, edge_weight: float) -> pd.DataFrame:
    rows = []
    for i, pid in enumerate(raw["patient_ids"]):
        dh = raw["delta_h"][i]
        da = raw["delta_alpha"][i]
        node_contrib = float(np.linalg.norm(dh, axis=1).mean())
        edge_contrib = float(np.abs(da).mean())
        transition = node_weight * node_contrib + edge_weight * edge_contrib

        row = {"patient_id": pid}
        for j, col in enumerate(exposure_cols):
            row[col] = float(raw["exposures"][i, j])
        row["NodeContribution"] = node_contrib
        row["EdgeContribution"] = edge_contrib
        row["TransitionScore"] = transition
        rows.append(row)

    return pd.DataFrame(rows).sort_values("TransitionScore", ascending=False).reset_index(drop=True)


def make_base_node_table(raw: Dict, node_labels: Dict[int, str], exposure_cols: List[str]) -> Tuple[pd.DataFrame, np.ndarray]:
    """NodeDelta is mean ||delta_h|| across patients. No smoking-correlation score is computed."""
    P = len(raw["patient_ids"])
    N = raw["n_nodes"]

    node_mag = np.zeros((P, N), dtype=np.float32)
    for i in range(P):
        node_mag[i] = np.linalg.norm(raw["delta_h"][i], axis=1)

    node_delta = node_mag.mean(axis=0)

    # Direction is only descriptive: larger node transition in high vs low primary exposure.
    # It is not used in the score.
    primary_exp = raw["exposures"][:, 0]
    high_mask = primary_exp >= np.median(primary_exp)
    low_mask = ~high_mask
    high_mean = node_mag[high_mask].mean(axis=0) if high_mask.any() else np.zeros(N)
    low_mean = node_mag[low_mask].mean(axis=0) if low_mask.any() else np.zeros(N)
    direction = np.sign(high_mean - low_mean).astype(int)

    rows = []
    for n in range(N):
        label = node_labels.get(n, f"node_{n}")
        node_type = node_type_from_index(n, raw["n_genes"])
        rows.append({
            "node_idx": int(n),
            "node_label": label,
            "node_type": node_type,
            "gene": clean_gene_label(label) if node_type == "gene" else "",
            "NodeDelta": float(node_delta[n]),
            "Direction": int(direction[n]),
        })

    return pd.DataFrame(rows), node_mag


def make_base_edge_table(raw: Dict, node_labels: Dict[int, str]) -> Tuple[pd.DataFrame, np.ndarray]:
    """EdgeDelta is mean |delta_alpha| across patients."""
    P = len(raw["patient_ids"])
    E = raw["edge_index"].shape[1]
    ei = raw["edge_index"]
    et = raw["edge_type"]

    da_abs = np.zeros((P, E), dtype=np.float32)
    da_raw = np.zeros((P, E), dtype=np.float32)
    alpha_n = np.zeros((P, E), dtype=np.float32)
    alpha_t = np.zeros((P, E), dtype=np.float32)

    for i in range(P):
        da_raw[i] = raw["delta_alpha"][i]
        da_abs[i] = np.abs(raw["delta_alpha"][i])
        alpha_n[i] = raw["alpha_n"][i]
        alpha_t[i] = raw["alpha_t"][i]

    edge_delta = da_abs.mean(axis=0)
    mean_delta_alpha = da_raw.mean(axis=0)
    direction = np.sign(mean_delta_alpha).astype(int)

    src = ei[0].astype(int)
    dst = ei[1].astype(int)

    rows = []
    for e in range(E):
        src_label = node_labels.get(int(src[e]), f"node_{int(src[e])}")
        dst_label = node_labels.get(int(dst[e]), f"node_{int(dst[e])}")
        rows.append({
            "edge_idx": int(e),
            "src_idx": int(src[e]),
            "src_label": src_label,
            "src_gene": clean_gene_label(src_label) if src_label.startswith("GENE:") else "",
            "dst_idx": int(dst[e]),
            "dst_label": dst_label,
            "dst_gene": clean_gene_label(dst_label) if dst_label.startswith("GENE:") else "",
            "edge_type_name": edge_type_name(et[e]),
            "EdgeDelta": float(edge_delta[e]),
            "BaseRewiringScore": float(edge_delta[e]),
            "mean_abs_delta_alpha": float(edge_delta[e]),
            "mean_delta_alpha": float(mean_delta_alpha[e]),
            "alpha_n_mean": float(alpha_n[:, e].mean()),
            "alpha_t_mean": float(alpha_t[:, e].mean()),
            "Direction": int(direction[e]),
        })

    edge_df = pd.DataFrame(rows).sort_values("EdgeDelta", ascending=False).reset_index(drop=True)
    return edge_df, da_abs


def add_network_weighted_node_score(
    node_df: pd.DataFrame,
    edge_df: pd.DataFrame,
    n_nodes: int,
    top_edge_fraction: float,
    node_delta_weight: float,
    weighted_degree_weight: float,
    smoking_assoc_weight: float,
    node_smoking_assoc: np.ndarray | None = None,
    smoking_assoc_col: str = "none",
    clustering_weight: float = 0.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    """Build high-rewiring network and compute final priority node scores."""
    if not 0.0 < top_edge_fraction < 1.0:
        raise ValueError("top_edge_fraction must be between 0 and 1")

    edge_cutoff = float(edge_df["EdgeDelta"].quantile(1.0 - top_edge_fraction))
    top_edges = edge_df[edge_df["EdgeDelta"] >= edge_cutoff].copy()
    top_edges["IsTopRewiredEdge"] = True
    top_edges["TopRewiredEdgeFraction"] = float(top_edge_fraction)
    top_edges["TopRewiredEdgeThreshold"] = edge_cutoff

    degree = np.zeros(n_nodes, dtype=int)
    weighted_degree = np.zeros(n_nodes, dtype=float)
    pairs: List[Tuple[int, int]] = []

    for _, row in top_edges.iterrows():
        s = int(row["src_idx"])
        d = int(row["dst_idx"])
        w = float(row["EdgeDelta"])
        if s == d:
            continue
        degree[s] += 1
        degree[d] += 1
        weighted_degree[s] += w
        weighted_degree[d] += w
        pairs.append((s, d))

    clustering = local_clustering_coefficients(n_nodes, pairs)

    scored = node_df.copy()
    scored["RewiredDegree"] = scored["node_idx"].astype(int).map(lambda x: int(degree[x]))
    scored["WeightedRewiredDegree"] = scored["node_idx"].astype(int).map(lambda x: float(weighted_degree[x]))
    scored["LocalClusteringCoefficient"] = scored["node_idx"].astype(int).map(lambda x: float(clustering[x]))

    # Diagnostic ranks are exported, but the final score below does NOT use ranks.
    scored["NodeDeltaRank"] = rank_norm(scored["NodeDelta"].values)
    scored["WeightedRewiredDegreeRank"] = rank_norm(scored["WeightedRewiredDegree"].values)
    scored["LocalClusteringRank"] = rank_norm(scored["LocalClusteringCoefficient"].values)

    # Final score used for TCGA pathway prioritization and simulation validation.
    # This uses magnitude-preserving min-max scaled components, not ranks.
    scored["NodeDeltaScaled"] = minmax01(scored["NodeDelta"].values)
    scored["WeightedRewiredDegreeScaled"] = minmax01(scored["WeightedRewiredDegree"].values)
    scored["LocalClusteringScaled"] = minmax01(scored["LocalClusteringCoefficient"].values)

    if node_smoking_assoc is None:
        node_smoking_assoc = np.zeros(len(scored), dtype=float)
    node_smoking_assoc = np.asarray(node_smoking_assoc, dtype=float)
    if len(node_smoking_assoc) != len(scored):
        raise ValueError(
            "node_smoking_assoc length does not match node table length: "
            f"{len(node_smoking_assoc)} vs {len(scored)}"
        )

    scored["NodeSmokingAssociationScore"] = node_smoking_assoc
    scored["NodeSmokingAssociationScoreScaled"] = minmax01(scored["NodeSmokingAssociationScore"].values)
    scored["NodeSmokingAssociationExposureCol"] = str(smoking_assoc_col)

    scored["FinalNodePriorityScore"] = (
        node_delta_weight * scored["NodeDeltaScaled"]
        + weighted_degree_weight * scored["WeightedRewiredDegreeScaled"]
        + smoking_assoc_weight * scored["NodeSmokingAssociationScoreScaled"]
        + clustering_weight * scored["LocalClusteringScaled"]
    )

    # Keep the variable name used by Step 11 / Step 14 / pathway scripts.
    scored["NetworkWeightedNodeScore"] = scored["FinalNodePriorityScore"]

    scored["StrictNetworkWeightedNodeScore"] = (
        scored["NodeDeltaScaled"] * np.sqrt(scored["WeightedRewiredDegreeScaled"] + 1e-12)
    )

    # Backward-compatible aliases for older plotting/evaluation scripts.
    scored["CancerNodeScore"] = scored["NodeDelta"]
    scored["CombinedNodeScore"] = scored["NetworkWeightedNodeScore"]

    # Top 10% flags for downstream pathway input.
    gene_mask = scored["node_type"].astype(str) == "gene"
    gene_scores = scored.loc[gene_mask, "NetworkWeightedNodeScore"]
    top_node_cutoff = float(gene_scores.quantile(0.90)) if len(gene_scores) else np.nan
    scored["IsTopNetworkWeightedGene10pct"] = False
    if np.isfinite(top_node_cutoff):
        scored.loc[gene_mask, "IsTopNetworkWeightedGene10pct"] = (
            scored.loc[gene_mask, "NetworkWeightedNodeScore"] >= top_node_cutoff
        )

    scored = scored.sort_values("NetworkWeightedNodeScore", ascending=False).reset_index(drop=True)
    return scored, top_edges.sort_values("EdgeDelta", ascending=False).reset_index(drop=True), edge_cutoff



def compute_signed_node_exposure_score(
    node_mag: np.ndarray,
    exposures: np.ndarray,
    exposure_cols: List[str],
    selected_exposure_col: str,
) -> np.ndarray:
    """Signed Spearman(node transition magnitude, selected exposure)."""
    if selected_exposure_col not in [str(c) for c in exposure_cols]:
        return np.zeros(node_mag.shape[1], dtype=float)
    exposure_idx = [str(c) for c in exposure_cols].index(selected_exposure_col)
    y = np.asarray(exposures[:, exposure_idx], dtype=float)
    scores = np.zeros(node_mag.shape[1], dtype=float)
    for j in range(node_mag.shape[1]):
        x = np.asarray(node_mag[:, j], dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 4 or np.nanstd(x[ok]) < 1e-12 or np.nanstd(y[ok]) < 1e-12:
            scores[j] = 0.0
            continue
        r, _ = spearmanr(x[ok], y[ok])
        scores[j] = float(r) if np.isfinite(r) else 0.0
    return scores


def apply_old_rank_based_edge_scores(
    edge_df: pd.DataFrame,
    raw: Dict,
    node_df: pd.DataFrame,
    exposure_cols: List[str],
    edge_scoring_mode: str = "exposure_high",
) -> pd.DataFrame:
    """
    Restore the older Step 05 edge scoring while keeping the current node scores.

    Node score columns remain unchanged. Edge final ranking uses rank-normalized
    tumor-normal rewiring, exposure association, exposure contrast, and endpoint
    support, exactly in the spirit of the older RewiringPotentialEdgeScore logic.

    Final edge ranking column:
      RewiringPotentialEdgeScore

    Compatibility:
      CombinedEdgeScore = RewiringPotentialEdgeScore
      EdgeDelta = RewiringPotentialEdgeScore   # so R scripts using EdgeDelta work
      RawEdgeDelta = original mean |delta_alpha|
    """
    P = len(raw["patient_ids"])
    E = raw["edge_index"].shape[1]
    if len(edge_df) != E:
        raise ValueError(f"edge_df length {len(edge_df)} != number of raw edges {E}")

    da_abs = np.zeros((P, E), dtype=np.float32)
    da_raw = np.zeros((P, E), dtype=np.float32)
    for i in range(P):
        da_raw[i] = raw["delta_alpha"][i]
        da_abs[i] = np.abs(raw["delta_alpha"][i])

    base_rewiring = da_abs.mean(axis=0)
    primary_exp = raw["exposures"][:, 0]
    high_mask = primary_exp >= np.median(primary_exp)
    low_mask = ~high_mask
    contrast = (
        da_abs[high_mask].mean(axis=0) if high_mask.any() else np.zeros(E)
    ) - (
        da_abs[low_mask].mean(axis=0) if low_mask.any() else np.zeros(E)
    )

    exp_scores = {}
    for j, col in enumerate(exposure_cols):
        exp = raw["exposures"][:, j]
        exp_scores[col] = np.array([
            spearmanr(da_abs[:, e], exp).statistic
            if (np.isfinite(da_abs[:, e]).sum() >= 4 and np.std(da_abs[:, e]) >= 1e-12 and np.std(exp) >= 1e-12)
            else 0.0
            for e in range(E)
        ], dtype=float)
        exp_scores[col] = np.nan_to_num(exp_scores[col], nan=0.0, posinf=0.0, neginf=0.0)

    mean_exp_score = np.stack(list(exp_scores.values()), axis=0).mean(axis=0) if exp_scores else np.zeros(E)

    rn_base = rank_norm(base_rewiring)
    rn_exp_pos = rank_norm(np.clip(mean_exp_score, 0.0, None))
    rn_contrast_pos = rank_norm(np.clip(contrast, 0.0, None))
    exposure_edge_signal = 0.5 * rn_exp_pos + 0.5 * rn_contrast_pos

    # Endpoint support from the current node table.
    n_nodes = raw["n_nodes"]
    node_idxed = node_df.set_index("node_idx", drop=False)
    cancer_node = np.zeros(n_nodes, dtype=float)
    exposure_node = np.zeros(n_nodes, dtype=float)

    for n in range(n_nodes):
        if n not in node_idxed.index:
            continue
        row = node_idxed.loc[n]
        cancer_node[n] = float(row.get("CancerNodeScore", row.get("NodeDelta", 0.0)))
        exp_cols_node = [c for c in node_df.columns if c.startswith("ExposureNodeScore_")]
        if exp_cols_node:
            exposure_node[n] = float(np.nanmean([row.get(c, 0.0) for c in exp_cols_node]))
        else:
            # Fallback to unsigned association if the signed column is unavailable.
            exposure_node[n] = float(row.get("NodeSmokingAssociationScore", 0.0))

    rn_cancer_node = rank_norm(cancer_node)
    rn_exposure_node = rank_norm(np.clip(exposure_node, 0.0, None))

    src = raw["edge_index"][0].astype(int)
    dst = raw["edge_index"][1].astype(int)

    exposure_endpoint_strict = np.sqrt(
        np.clip(rn_exposure_node[src], 0.0, None) * np.clip(rn_exposure_node[dst], 0.0, None)
    )
    exposure_endpoint_anchor = np.maximum(rn_exposure_node[src], rn_exposure_node[dst])
    cancer_endpoint_support = np.sqrt(
        np.clip(rn_cancer_node[src], 0.0, None) * np.clip(rn_cancer_node[dst], 0.0, None)
    )

    if edge_scoring_mode == "balanced":
        cancer_w, exposure_w = 0.5, 0.5
    elif edge_scoring_mode == "exposure_high":
        cancer_w, exposure_w = 0.3, 0.7
    elif edge_scoring_mode == "exposure_only":
        cancer_w, exposure_w = 0.0, 1.0
    else:
        raise ValueError("edge_scoring_mode must be one of: balanced, exposure_high, exposure_only")

    mixed_node_support = cancer_w * rn_cancer_node + exposure_w * rn_exposure_node
    mixed_endpoint_support = np.sqrt(
        np.clip(mixed_node_support[src], 0.0, None) * np.clip(mixed_node_support[dst], 0.0, None)
    )

    cancer_rewiring_score = rn_base * cancer_endpoint_support
    exposure_dependent_score = exposure_edge_signal * exposure_endpoint_strict
    exposure_anchored_score = exposure_edge_signal * exposure_endpoint_anchor
    mixed_cancer_exposure_score = (0.5 * rn_base + 0.5 * exposure_edge_signal) * mixed_endpoint_support

    if edge_scoring_mode == "balanced":
        final_score = mixed_cancer_exposure_score
        final_score_source = "MixedCancerExposureRewiringScore"
    elif edge_scoring_mode == "exposure_high":
        final_score = exposure_anchored_score
        final_score_source = "ExposureAnchoredRewiringEdgeScore"
    else:
        final_score = exposure_dependent_score
        final_score_source = "ExposureDependentRewiringScore"

    out = edge_df.copy()
    out["RawEdgeDelta"] = base_rewiring
    out["BaseRewiringScore"] = base_rewiring
    out["mean_abs_delta_alpha"] = base_rewiring
    out["ExposureContrastScore"] = contrast
    out["ExposureEdgeSignal"] = exposure_edge_signal
    out["CancerEndpointSupport"] = cancer_endpoint_support
    out["ExposureEndpointSupportStrict"] = exposure_endpoint_strict
    out["ExposureEndpointSupportAnchor"] = exposure_endpoint_anchor
    out["MixedEndpointSupport"] = mixed_endpoint_support
    out["CancerRewiringPotentialEdgeScore"] = cancer_rewiring_score
    out["ExposureDependentRewiringScore"] = exposure_dependent_score
    out["ExposureAnchoredRewiringEdgeScore"] = exposure_anchored_score
    out["MixedCancerExposureRewiringScore"] = mixed_cancer_exposure_score

    for col in exposure_cols:
        out[f"ExposureRewiringScore_{col}"] = exp_scores[col]

    out["RewiringPotentialEdgeScore"] = final_score
    out["CombinedEdgeScore"] = out["RewiringPotentialEdgeScore"]
    # Important compatibility choice: R scripts that rank by EdgeDelta now get
    # the same final edge score. Raw edge delta is kept in RawEdgeDelta.
    out["EdgeDelta"] = out["RewiringPotentialEdgeScore"]
    out["EdgeScoringMode"] = str(edge_scoring_mode)
    out["RewiringPotentialEdgeScoreSource"] = final_score_source

    out = out.sort_values("RewiringPotentialEdgeScore", ascending=False).reset_index(drop=True)

    print("\nFinal edge score used for ranking:")
    print(f"  edge_scoring_mode = {edge_scoring_mode}")
    print(f"  RewiringPotentialEdgeScore source = {final_score_source}")
    print(
        "  score min/median/max = "
        f"{out['RewiringPotentialEdgeScore'].min():.4f} / "
        f"{out['RewiringPotentialEdgeScore'].median():.4f} / "
        f"{out['RewiringPotentialEdgeScore'].max():.4f}"
    )

    return out

def make_pathway_ready_tables(node_df: pd.DataFrame, top_edges_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Create simple gene lists for downstream node- and edge-based pathway analysis."""
    node_pathway_input = node_df[
        (node_df["node_type"] == "gene") & (node_df["IsTopNetworkWeightedGene10pct"])
    ].copy()

    node_pathway_input = node_pathway_input[[
        "gene",
        "node_idx",
        "node_label",
        "NodeDelta",
        "RewiredDegree",
        "WeightedRewiredDegree",
        "LocalClusteringCoefficient",
        "NodeSmokingAssociationScore",
        "FinalNodePriorityScore",
        "NetworkWeightedNodeScore",
    ]].sort_values("NetworkWeightedNodeScore", ascending=False).reset_index(drop=True)

    edge_genes = []
    for _, row in top_edges_df.iterrows():
        if str(row["src_label"]).startswith("GENE:"):
            edge_genes.append(row["src_gene"])
        if str(row["dst_label"]).startswith("GENE:"):
            edge_genes.append(row["dst_gene"])

    edge_pathway_input = pd.DataFrame({"gene": sorted(set(g for g in edge_genes if str(g) != ""))})
    edge_pathway_input = edge_pathway_input.merge(
        node_df[[
            "gene",
            "node_idx",
            "node_label",
            "NodeDelta",
            "RewiredDegree",
            "WeightedRewiredDegree",
            "LocalClusteringCoefficient",
            "NodeSmokingAssociationScore",
            "FinalNodePriorityScore",
            "NetworkWeightedNodeScore",
        ]],
        on="gene",
        how="left",
    ).sort_values("NetworkWeightedNodeScore", ascending=False).reset_index(drop=True)

    return node_pathway_input, edge_pathway_input


def make_patient_rewiring_burden(
    raw: Dict,
    edge_abs_matrix: np.ndarray,
    top_edges_df: pd.DataFrame,
    exposure_cols: List[str],
) -> pd.DataFrame:
    selected = top_edges_df["edge_idx"].astype(int).values
    weights = top_edges_df["EdgeDelta"].astype(float).values

    rows = []
    for i, pid in enumerate(raw["patient_ids"]):
        row = {"patient_id": pid}
        for j, col in enumerate(exposure_cols):
            row[col] = float(raw["exposures"][i, j])

        if len(selected) == 0:
            raw_burden = 0.0
            weighted_burden = 0.0
        else:
            vals = edge_abs_matrix[i, selected]
            raw_burden = float(np.mean(vals))
            weighted_burden = float(np.sum(vals * weights) / max(np.sum(weights), 1e-12))

        row["SelectedTopRewiredEdgeCount"] = int(len(selected))
        row["PatientRewiringBurden"] = raw_burden
        row["WeightedPatientRewiringBurden"] = weighted_burden
        rows.append(row)

    return pd.DataFrame(rows).sort_values("WeightedPatientRewiringBurden", ascending=False).reset_index(drop=True)


def make_per_patient_gene_matrix(raw: Dict, node_labels: Dict[int, str], exposure_cols: List[str], node_mag: np.ndarray) -> pd.DataFrame:
    col_labels = [node_labels.get(n, f"node_{n}") for n in range(raw["n_nodes"])]
    df = pd.DataFrame(node_mag, columns=col_labels)
    df.insert(0, "patient_id", raw["patient_ids"])
    for j, col in enumerate(exposure_cols):
        df.insert(j + 1, col, raw["exposures"][:, j])
    return df


def make_per_patient_edge_matrix(raw: Dict, exposure_cols: List[str], edge_abs_matrix: np.ndarray) -> pd.DataFrame:
    df = pd.DataFrame(edge_abs_matrix, columns=[f"edge_{e}" for e in range(edge_abs_matrix.shape[1])])
    df.insert(0, "patient_id", raw["patient_ids"])
    for j, col in enumerate(exposure_cols):
        df.insert(j + 1, col, raw["exposures"][:, j])
    return df


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    ap = argparse.ArgumentParser(description="Clean ExposureGAT result extraction")
    ap.add_argument("--data_dir", type=str, default="GNN_Input_Data")
    ap.add_argument("--model_path", type=str, default="GNN_Output/best_model.pth")
    ap.add_argument("--config", type=str, default="GNN_Output/run_config.json")
    ap.add_argument("--out_dir", type=str, default="GNN_Output/results")
    ap.add_argument("--node_weight", type=float, default=1.0)
    ap.add_argument("--edge_weight", type=float, default=1.0)
    ap.add_argument("--top_nodes", type=int, default=200)
    ap.add_argument("--top_edges", type=int, default=500)
    ap.add_argument("--top_edge_fraction", type=float, default=0.10)
    # Final node priority score. Defaults:
    #   0.70 * scaled(NodeDelta)
    # + 0.20 * scaled(WeightedRewiredDegree)
    # + 0.10 * scaled(NodeSmokingAssociationScore)
    # Clustering is kept as an optional extra term for sensitivity analysis.
    ap.add_argument("--node_delta_weight", type=float, default=0.70)
    ap.add_argument("--weighted_degree_weight", type=float, default=0.20)
    ap.add_argument("--smoking_assoc_weight", type=float, default=0.10)
    ap.add_argument("--clustering_weight", type=float, default=0.00)
    ap.add_argument("--smoking_assoc_col", type=str, default=None,
                    help="Exposure column to use for NodeSmokingAssociationScore. Defaults to smoking_intensity, pack_years, or first exposure column.")

    # Backward-compatible ignored arguments from the old messy Step 05.
    # Keeping them prevents old run commands from crashing, but they no longer do anything here.
    ap.add_argument("--gmt_files", type=str, nargs="+", default=None)
    ap.add_argument("--min_coverage_frac", type=float, default=0.10)
    ap.add_argument("--min_genes_in_model", type=int, default=10)
    ap.add_argument("--top_nodes_assoc", type=int, default=30)
    ap.add_argument("--top_edges_per_node", type=int, default=5)
    ap.add_argument("--rewired_edge_top_fraction", type=float, default=None)
    ap.add_argument("--edge_top_percent", type=float, default=None)
    ap.add_argument(
        "--edge_scoring_mode",
        type=str,
        default="exposure_high",
        choices=["balanced", "exposure_high", "exposure_only"],
        help="Old rank-based final edge score: balanced=mixed cancer+exposure; exposure_high=exposure-anchored; exposure_only=strict exposure-dependent.",
    )

    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # If the old argument is supplied, let it control the top-edge network fraction.
    if args.edge_top_percent is not None:
        args.top_edge_fraction = args.edge_top_percent
    elif args.rewired_edge_top_fraction is not None:
        args.top_edge_fraction = args.rewired_edge_top_fraction

    weight_sum = (
        args.node_delta_weight
        + args.weighted_degree_weight
        + args.smoking_assoc_weight
        + args.clustering_weight
    )
    if not np.isclose(weight_sum, 1.0):
        print(f"WARNING: node score weights sum to {weight_sum:.3f}, not 1.0. Proceeding anyway.")

    if args.gmt_files:
        print("NOTE: --gmt_files was provided, but pathway analysis has been removed from Step 05.")
        print("      Use Step 06 for pathway analysis using Step 05 exported pathway input files.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    with open(args.config) as f:
        cfg = json.load(f)
    print(f"Config loaded from {args.config}")

    dataset = LungMultiOmicsDataset(
        data_dir=args.data_dir,
        exposure_cols=tuple(cfg.get("exposure_cols", ["smoking_intensity"])),
        meth_min_scale=cfg.get("meth_min_scale", 0.005),
        cnv_min_scale=cfg.get("cnv_min_scale", 0.05),
        verbose=True,
    )
    dataset.fit_scalers(train_patient_ids=dataset.patient_ids)

    if cfg.get("shuffle_exposure", False):
        print("\n" + "!" * 60)
        print("WARNING: config is from a SHUFFLED EXPOSURE control run.")
        print("Do NOT use these results for biological interpretation.")
        print("!" * 60 + "\n")

    if cfg.get("no_exposure_model", False):
        print("\n" + "!" * 60)
        print("WARNING: config is from a NO-EXPOSURE control run.")
        print("Scores reflect graph-only transition, not exposure-gated rewiring.")
        print("!" * 60 + "\n")

    n_nodes_cfg = cfg.get("n_nodes", dataset.n_nodes)
    if n_nodes_cfg != dataset.n_nodes:
        print(f"WARNING: config n_nodes={n_nodes_cfg} differs from dataset n_nodes={dataset.n_nodes}; using dataset value.")
        n_nodes_cfg = dataset.n_nodes

    model = ExposureGAT(
        cpg_in=1,
        gene_in=3,
        exp_dim=len(dataset.exposure_cols),
        hidden_dim=cfg.get("hidden_dim", 64),
        out_dim=cfg.get("out_dim", 32),
        heads=cfg.get("heads", 4),
        dropout=0.0,
        tanh_scale=cfg.get("tanh_scale", 0.5),
        bias_hidden=cfg.get("bias_hidden", 32),
        nudge_scale=cfg.get("nudge_scale", 0.1),
        n_nodes=n_nodes_cfg,
        edge_index_for_degree=dataset.edge_index,
    )

    try:
        state = torch.load(args.model_path, map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(args.model_path, map_location=device)
    model.load_state_dict(state)
    model.eval()
    print(f"Model loaded from {args.model_path}")

    node_labels = dataset.get_node_labels()
    exposure_cols = list(dataset.exposure_cols)

    print("\nRunning inference...")
    raw = run_inference(model, dataset, device)

    print("\nBuilding patient table...")
    patient_df = make_patient_table(raw, exposure_cols, args.node_weight, args.edge_weight)

    print("Building base node table...")
    node_df, node_mag = make_base_node_table(raw, node_labels, exposure_cols)

    print("Building base edge table...")
    edge_df, edge_abs_matrix = make_base_edge_table(raw, node_labels)

    print("Computing node smoking/exposure association score...")
    node_smoking_assoc, smoking_assoc_col = compute_node_smoking_association(
        node_mag=node_mag,
        exposures=raw["exposures"],
        exposure_cols=exposure_cols,
        preferred_exposure_col=args.smoking_assoc_col,
    )
    print(f"Using exposure column for NodeSmokingAssociationScore: {smoking_assoc_col}")

    # Add signed exposure association for compatibility with the older edge scoring logic.
    # NodeSmokingAssociationScore remains absolute for the current node score.
    signed_node_exposure = compute_signed_node_exposure_score(
        node_mag=node_mag,
        exposures=raw["exposures"],
        exposure_cols=exposure_cols,
        selected_exposure_col=smoking_assoc_col,
    )
    node_df[f"ExposureNodeScore_{smoking_assoc_col}"] = signed_node_exposure

    print("Building high-rewiring network and final priority node score...")
    node_df, top_edges_df, edge_cutoff = add_network_weighted_node_score(
        node_df=node_df,
        edge_df=edge_df,
        n_nodes=raw["n_nodes"],
        top_edge_fraction=args.top_edge_fraction,
        node_delta_weight=args.node_delta_weight,
        weighted_degree_weight=args.weighted_degree_weight,
        smoking_assoc_weight=args.smoking_assoc_weight,
        node_smoking_assoc=node_smoking_assoc,
        smoking_assoc_col=smoking_assoc_col,
        clustering_weight=args.clustering_weight,
    )

    # Node scores above are intentionally unchanged from the current preferred version.
    # Now restore the older rank-based edge scoring for edge ranking/plots.
    print("Applying old rank-based RewiringPotentialEdgeScore edge logic...")
    edge_df = apply_old_rank_based_edge_scores(
        edge_df=edge_df,
        raw=raw,
        node_df=node_df,
        exposure_cols=exposure_cols,
        edge_scoring_mode=args.edge_scoring_mode,
    )

    # Recompute top edges using the restored final edge score.
    edge_cutoff = float(edge_df["RewiringPotentialEdgeScore"].quantile(1.0 - args.top_edge_fraction))
    edge_df["IsTopRewiredEdge"] = edge_df["RewiringPotentialEdgeScore"] >= edge_cutoff
    edge_df["TopRewiredEdgeFraction"] = float(args.top_edge_fraction)
    edge_df["TopRewiredEdgeThreshold"] = float(edge_cutoff)
    top_edges_df = edge_df[edge_df["IsTopRewiredEdge"]].copy().sort_values(
        "RewiringPotentialEdgeScore", ascending=False
    ).reset_index(drop=True)
    edge_df = edge_df.sort_values("RewiringPotentialEdgeScore", ascending=False).reset_index(drop=True)

    print("Building pathway-ready gene lists...")
    node_pathway_input, edge_pathway_input = make_pathway_ready_tables(node_df, top_edges_df)

    print("Building patient rewiring burden table...")
    patient_burden_df = make_patient_rewiring_burden(raw, edge_abs_matrix, top_edges_df, exposure_cols)

    print("Building per-patient matrices for plotting...")
    per_patient_gene_df = make_per_patient_gene_matrix(raw, node_labels, exposure_cols, node_mag)
    per_patient_edge_df = make_per_patient_edge_matrix(raw, exposure_cols, edge_abs_matrix)

    # Compact reporting outputs.
    node_report_cols = [
        "node_idx",
        "node_label",
        "node_type",
        "gene",
        "NodeDelta",
        "Direction",
        "RewiredDegree",
        "WeightedRewiredDegree",
        "LocalClusteringCoefficient",
        "NodeDeltaScaled",
        "WeightedRewiredDegreeScaled",
        "LocalClusteringScaled",
        "NodeSmokingAssociationScore",
        "NodeSmokingAssociationScoreScaled",
        "NodeSmokingAssociationExposureCol",
        f"ExposureNodeScore_{smoking_assoc_col}",
        "FinalNodePriorityScore",
        "NodeDeltaRank",
        "WeightedRewiredDegreeRank",
        "LocalClusteringRank",
        "NetworkWeightedNodeScore",
        "StrictNetworkWeightedNodeScore",
        "CancerNodeScore",
        "CombinedNodeScore",
        "IsTopNetworkWeightedGene10pct",
    ]

    edge_report_cols = [
        "edge_idx",
        "src_idx",
        "src_label",
        "src_gene",
        "dst_idx",
        "dst_label",
        "dst_gene",
        "edge_type_name",
        "EdgeDelta",
        "RawEdgeDelta",
        "BaseRewiringScore",
        "ExposureEdgeSignal",
        "ExposureContrastScore",
        "CancerRewiringPotentialEdgeScore",
        "ExposureDependentRewiringScore",
        "ExposureAnchoredRewiringEdgeScore",
        "MixedCancerExposureRewiringScore",
        "RewiringPotentialEdgeScore",
        "CombinedEdgeScore",
        "mean_abs_delta_alpha",
        "mean_delta_alpha",
        "alpha_n_mean",
        "alpha_t_mean",
        "Direction",
        "IsTopRewiredEdge",
        "TopRewiredEdgeFraction",
        "TopRewiredEdgeThreshold",
    ]

    # top_edges_df was created before backward-compatible edge aliases were added
    # to edge_df. Rebuild the top-edge report from the full edge_df so all
    # report columns are guaranteed to exist.
    top_edges_report = edge_df.loc[edge_df["IsTopRewiredEdge"].astype(bool), edge_report_cols].copy()

    summary_df = pd.DataFrame({
        "parameter": [
            "n_patients",
            "n_nodes",
            "n_gene_nodes",
            "n_edges",
            "top_edge_fraction",
            "top_edge_delta_threshold",
            "n_top_rewired_edges",
            "node_delta_weight",
            "weighted_degree_weight",
            "smoking_assoc_weight",
            "clustering_weight",
            "smoking_assoc_col",
            "n_node_pathway_genes",
            "n_edge_pathway_genes",
            "node_score_formula",
            "edge_score_formula",
            "edge_scoring_mode",
        ],
        "value": [
            len(raw["patient_ids"]),
            raw["n_nodes"],
            raw["n_genes"],
            raw["edge_index"].shape[1],
            args.top_edge_fraction,
            edge_cutoff,
            len(top_edges_df),
            args.node_delta_weight,
            args.weighted_degree_weight,
            args.smoking_assoc_weight,
            args.clustering_weight,
            smoking_assoc_col,
            len(node_pathway_input),
            len(edge_pathway_input),
            "NetworkWeightedNodeScore = FinalNodePriorityScore = w1*scaled(NodeDelta) + w2*scaled(WeightedRewiredDegree) + w3*scaled(NodeSmokingAssociationScore) + optional w4*scaled(LocalClusteringCoefficient); scaled = min-max 0-1, not rank",
            "RewiringPotentialEdgeScore restored from old rank-based formula; EdgeDelta is set equal to RewiringPotentialEdgeScore for R compatibility; RawEdgeDelta stores mean(abs(delta_alpha))",
            args.edge_scoring_mode,
        ],
    })

    print("\nSaving clean Step 05 outputs...")
    write_csv(patient_df, args.out_dir, "results_patients.csv")
    write_csv(node_df[node_report_cols], args.out_dir, "results_nodes.csv")
    write_csv(node_df[node_report_cols].head(args.top_nodes), args.out_dir, f"results_nodes_top{args.top_nodes}.csv")
    write_csv(edge_df[edge_report_cols], args.out_dir, "results_edges.csv")
    write_csv(edge_df[edge_report_cols].head(args.top_edges), args.out_dir, f"results_edges_top{args.top_edges}.csv")
    write_csv(top_edges_report, args.out_dir, "results_top_rewired_edges.csv")
    write_csv(node_pathway_input, args.out_dir, "pathway_input_node_top10pct_network_weighted_genes.csv")
    write_csv(edge_pathway_input, args.out_dir, "pathway_input_edge_genes_touched_by_top_rewired_edges.csv")
    write_csv(patient_burden_df, args.out_dir, "results_patient_rewiring_burden.csv")
    write_csv(per_patient_gene_df, args.out_dir, "results_per_patient_gene_scores.csv")
    write_csv(per_patient_edge_df, args.out_dir, "results_per_patient_edge_rewiring.csv")
    write_csv(summary_df, args.out_dir, "step05_clean_summary.csv")

    print("\n" + "=" * 70)
    print("STEP 05 CLEAN RESULTS SUMMARY")
    print("=" * 70)
    print(f"Patients processed: {len(raw['patient_ids'])}")
    print(f"Nodes processed: {raw['n_nodes']}  gene nodes: {raw['n_genes']}")
    print(f"Edges processed: {raw['edge_index'].shape[1]}")
    print(f"Top rewired-edge fraction: {args.top_edge_fraction:.1%}")
    print(f"EdgeDelta threshold: {edge_cutoff:.6f}")
    print("\nNode score:")
    print("  FinalNodePriorityScore = NetworkWeightedNodeScore =")
    print(f"    {args.node_delta_weight:.2f} * scaled(NodeDelta)")
    print(f"  + {args.weighted_degree_weight:.2f} * scaled(WeightedRewiredDegree)")
    print(f"  + {args.smoking_assoc_weight:.2f} * scaled(NodeSmokingAssociationScore)")
    print(f"  + {args.clustering_weight:.2f} * scaled(LocalClusteringCoefficient)")
    print(f"  NodeSmokingAssociationScore exposure column: {smoking_assoc_col}")
    print("  where scaled = min-max 0-1 magnitude scaling, not rank")
    print("\nTop 15 final-priority genes:")
    show_nodes = node_df[node_df["node_type"] == "gene"][[
        "node_label", "NodeDelta", "WeightedRewiredDegree", "NodeSmokingAssociationScore", "NetworkWeightedNodeScore"
    ]].head(15)
    print(show_nodes.to_string(index=False))
    print("\nTop 15 rewired edges by EdgeDelta:")
    show_edges = edge_df[[
        "src_label", "dst_label", "edge_type_name", "EdgeDelta", "Direction"
    ]].head(15)
    print(show_edges.to_string(index=False))
    print("=" * 70)
    print("\nPathway analysis is intentionally removed from Step 05.")
    print("Use these files in Step 06:")
    print("  pathway_input_node_top10pct_network_weighted_genes.csv")
    print("  pathway_input_edge_genes_touched_by_top_rewired_edges.csv")


if __name__ == "__main__":
    main()
