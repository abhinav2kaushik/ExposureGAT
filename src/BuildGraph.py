#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuildGraph.py

Dataset loader for ExposureGAT.

Expected files in data_dir (written by Step_01):
  X_rna_gene_features.csv   [samples x genes]  log1p RNA expression
  X_meth_gene_features.csv  [samples x genes]  mean promoter beta per gene
  X_cnv_gene_features.csv   [samples x genes]  sign*log2(|copy_number-2|+1)
  X_cpg_features.csv        [samples x cpgs]   raw beta per CpG node
  edges_ppi.csv             from_idx, to_idx, weight
  edges_cpg_to_gene.csv     from_idx, to_idx, weight, cpg_id, entrez_id
  map_gene_nodes.csv        entrez_id, node_idx, node_type
  map_cpg_nodes.csv         cpg_id,    node_idx, node_type
  metadata_samples.csv      sample_id, patient_id, tissue_type,
                            smoking_intensity, age, gender, ...

Node layout (fixed across all samples):
  indices 0 .. n_genes-1          -> gene nodes  (node_type = 1)
  indices n_genes .. n_nodes-1    -> CpG  nodes  (node_type = 0)

Node features packed into x tensor of shape [N, 3]:
  gene node : [RNA, mean_promoter_beta, CNV]   columns 0, 1, 2
  CpG  node : [beta, 0, 0]                     column 0 only
  NodeTypeProjection in Step_03 reads each type using only its own columns.

Edges:
  PPI (gene<->gene): undirected, both directions already in CSV from Step_01.
  CpG<->Gene: Step_01 writes CpG->Gene only. This loader automatically adds
  the reverse Gene->CpG direction so CpG nodes receive gene neighbourhood
  context during message passing, enabling methylation embeddings to reflect
  the expression/CNV state of the regulated gene.

Scaling (fit on training patients only, no leakage):
  Robust scaling: center = median, scale = max(IQR, min_scale_floor).
  Followed by tanh to suppress outliers.
  CpG beta values are NOT scaled (already in [0, 1]).
  CNV min_scale = 0.05 (log2 units from Step_01 transform).
  Meth min_scale = 0.01 (low floor lets real variation through).

Returns per __getitem__:
  (normal_graph, tumor_graph, exposure, patient_id)
  exposure : [exp_dim] — standardised pack_years (zero mean, unit std
             on training patients).  Raw values are NOT returned.
             Standardisation prevents tanh saturation in ExposureEdgeBias
             and ExposureConditioner MLPs for high pack-years patients.
             Parameters (_exp_mean, _exp_std) are fit in fit_scalers()
             on training patients only — no leakage into validation.
"""

import os
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data


class MultiOmicsDataset(Dataset):

    def __init__(
        self,
        data_dir:       str           = "GNN_Input_Data",
        exposure_cols:  Sequence[str] = ("smoking_intensity",),
        sample_id_col:  str           = "sample_id",
        patient_id_col: str           = "patient_id",
        tissue_col:     str           = "tissue_type",
        tumor_label:    Sequence[str] = ("TP",),
        normal_label:   Sequence[str] = ("NT",),
        rna_min_scale:  float         = 0.5,
        meth_min_scale: float         = 0.005,
        cnv_min_scale:  float         = 0.05,
        tanh_bound:     float         = 3.0,
        verbose:        bool          = True,
    ):
        self.data_dir       = data_dir
        self.exposure_cols  = list(exposure_cols)
        self.sample_id_col  = sample_id_col
        self.patient_id_col = patient_id_col
        self.tissue_col     = tissue_col
        self.tumor_label    = tuple(tumor_label)
        self.normal_label   = tuple(normal_label)
        self.rna_min_scale  = rna_min_scale
        self.meth_min_scale = meth_min_scale
        self.cnv_min_scale  = cnv_min_scale
        self.tanh_bound     = tanh_bound
        self.verbose        = verbose

        self._load_features()
        self._load_metadata()
        self._load_edges()
        self._build_node_type_tensor()

        self.scalers_fitted  = False
        self._rna_center     = self._rna_scale  = None
        self._meth_center    = self._meth_scale = None
        self._cnv_center     = self._cnv_scale  = None

        if self.verbose:
            self._print_summary()

    # ------------------------------------------------------------------
    # Internal loading
    # ------------------------------------------------------------------

    def _read_csv(self, fname: str) -> pd.DataFrame:
        df = pd.read_csv(os.path.join(self.data_dir, fname), index_col=0)
        df.index   = df.index.astype(str)
        df.columns = df.columns.astype(str)
        return df.astype(np.float32)

    def _load_features(self) -> None:
        self.rna  = self._read_csv("X_rna_gene_features.csv")
        self.meth = self._read_csv("X_meth_gene_features.csv")
        self.cnv  = self._read_csv("X_cnv_gene_features.csv")
        self.cpg  = self._read_csv("X_cpg_features.csv")

        for name, df in [("meth", self.meth), ("cnv", self.cnv), ("cpg", self.cpg)]:
            if not self.rna.index.equals(df.index):
                raise ValueError(
                    f"Sample row order mismatch between RNA and {name}. "
                    "All feature CSVs must have identical row order."
                )

        if not (self.rna.columns.equals(self.meth.columns) and
                self.rna.columns.equals(self.cnv.columns)):
            raise ValueError(
                "Gene column order mismatch across RNA / meth / CNV matrices."
            )

        self.genes   = list(self.rna.columns)
        self.cpgs    = list(self.cpg.columns)
        self.n_genes = len(self.genes)
        self.n_cpgs  = len(self.cpgs)
        self.n_nodes = self.n_genes + self.n_cpgs

    def _load_metadata(self) -> None:
        meta = pd.read_csv(os.path.join(self.data_dir, "metadata_samples.csv"))

        # Accept "sample_id" or "sample_barcode" as the sample column
        if self.sample_id_col not in meta.columns:
            for candidate in ("sample_id", "sample_barcode"):
                if candidate in meta.columns:
                    self.sample_id_col = candidate
                    break
            else:
                raise ValueError(
                    "metadata_samples.csv must contain 'sample_id' or 'sample_barcode'."
                )

        required = ([self.sample_id_col, self.patient_id_col, self.tissue_col]
                    + self.exposure_cols)
        missing = [c for c in required if c not in meta.columns]
        if missing:
            raise ValueError(f"metadata_samples.csv missing columns: {missing}")

        meta[self.sample_id_col]  = meta[self.sample_id_col].astype(str)
        meta[self.patient_id_col] = meta[self.patient_id_col].astype(str)

        def _classify(val):
            s = str(val).upper()
            if any(t.upper() in s for t in self.tumor_label):
                return "tumor"
            if any(t.upper() in s for t in self.normal_label):
                return "normal"
            return None

        meta["_tissue"] = meta[self.tissue_col].apply(_classify)

        bad = meta["_tissue"].isna()
        if bad.any():
            examples = meta.loc[bad, self.tissue_col].unique()[:5]
            raise ValueError(
                f"Could not classify tissue values as tumor/normal: {examples}"
            )

        # Keep only samples present in feature matrices
        feat_samples = set(self.rna.index)
        meta = meta[meta[self.sample_id_col].isin(feat_samples)].copy()

        # Drop patients with missing exposure values
        meta = meta.dropna(subset=self.exposure_cols).copy()

        # Keep only patients with both a normal and a tumor sample
        has_both = (
            meta.groupby(self.patient_id_col)["_tissue"]
            .apply(lambda x: {"normal", "tumor"}.issubset(set(x)))
        )
        valid_pids = has_both[has_both].index.tolist()
        meta = meta[meta[self.patient_id_col].isin(valid_pids)].copy()

        # One row per (patient, tissue) — first by barcode if duplicates
        meta = (
            meta
            .sort_values([self.patient_id_col, "_tissue", self.sample_id_col])
            .groupby([self.patient_id_col, "_tissue"], as_index=False)
            .first()
        )

        self.meta        = meta
        self.patient_ids = sorted(meta[self.patient_id_col].unique().tolist())

    def _load_edges(self) -> None:
        """
        Build shared edge_index [2, E], edge_attr [E, 1], edge_type [E].

        PPI edges: undirected, both directions in edges_ppi.csv.  edge_type = 0
        CpG-Gene : Step_01 writes CpG->Gene only; Gene->CpG added here. edge_type = 1

        edge_type is passed to Step_03 ExposureEdgeBias so the smoking bias
        MLP can learn separate modulation functions for PPI vs CpG-Gene edges.
        This is biologically motivated: smoking may differentially rewire
        epigenetic regulation (CpG-Gene) vs protein interaction (PPI) circuits.
        """
        ppi = pd.read_csv(os.path.join(self.data_dir, "edges_ppi.csv"))
        cpg = pd.read_csv(os.path.join(self.data_dir, "edges_cpg_to_gene.csv"))

        ppi_ei = ppi[["from_idx", "to_idx"]].to_numpy(dtype=np.int64)
        cpg_ei = cpg[["from_idx", "to_idx"]].to_numpy(dtype=np.int64)

        # Apply n_genes offset to CpG source indices if Step_01 didn't.
        # Guard: skip if CpG edge file is empty (e.g. simulation with no CpG nodes).
        if len(cpg_ei) > 0 and cpg_ei[:, 0].min() < self.n_genes:
            cpg_ei[:, 0] += self.n_genes

        cpg_ei_rev = cpg_ei[:, [1, 0]] if len(cpg_ei) > 0 else cpg_ei   # Gene->CpG reverse

        edge_index = np.vstack([ppi_ei, cpg_ei, cpg_ei_rev]) if len(cpg_ei) > 0 \
                     else ppi_ei
        n_edges    = len(edge_index)

        if edge_index.min() < 0 or edge_index.max() >= self.n_nodes:
            raise ValueError(
                f"Edge indices out of range [0, {self.n_nodes}). "
                f"Got min={edge_index.min()}, max={edge_index.max()}."
            )

        self.edge_index = torch.tensor(
            edge_index, dtype=torch.long
        ).t().contiguous()
        self.edge_attr  = torch.ones(n_edges, 1, dtype=torch.float32)

        # edge_type: 0=PPI, 1=CpG-Gene (both directions)
        n_ppi = len(ppi_ei)
        n_cpg = len(cpg_ei) * 2   # forward + reverse
        self.edge_type = torch.cat([
            torch.zeros(n_ppi, dtype=torch.long),
            torch.ones(n_cpg,  dtype=torch.long),
        ])   # [E]

        if self.verbose:
            print(f"  edges_ppi           : {n_ppi}")
            print(f"  edges_CpG<->Gene    : {len(cpg_ei)} each dir "
                  f"(={n_cpg} total)")
            print(f"  total edges         : {n_edges}")

    def _build_node_type_tensor(self) -> None:
        """node_type = 1 for gene nodes, 0 for CpG nodes."""
        types = torch.zeros(self.n_nodes, dtype=torch.long)
        types[:self.n_genes] = 1
        self.node_type = types

    # ------------------------------------------------------------------
    # Scaling — fit on training patients only
    # ------------------------------------------------------------------

    @staticmethod
    def _robust_scale(
        X: np.ndarray, min_scale: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        center = np.median(X, axis=0).astype(np.float32)
        q75    = np.percentile(X, 75, axis=0).astype(np.float32)
        q25    = np.percentile(X, 25, axis=0).astype(np.float32)
        scale  = np.maximum(q75 - q25, min_scale).astype(np.float32)
        return center, scale

    def fit_scalers(self, train_patient_ids: Iterable[str]) -> None:
        """
        Fit robust scalers on training patients only.
        Must be called after the train/val split and before __getitem__.
        """
        train_pids    = set(str(p) for p in train_patient_ids)
        train_meta    = self.meta[self.meta[self.patient_id_col].isin(train_pids)]
        train_samples = [s for s in train_meta[self.sample_id_col].astype(str).tolist()
                         if s in self.rna.index]

        if not train_samples:
            raise ValueError(
                "fit_scalers: no training samples matched feature matrix rows."
            )

        rna_tr  = self.rna.loc[train_samples].to_numpy(dtype=np.float32)
        meth_tr = self.meth.loc[train_samples].to_numpy(dtype=np.float32)
        cnv_tr  = np.nan_to_num(
            self.cnv.loc[train_samples].to_numpy(dtype=np.float32), nan=0.0
        )

        self._rna_center,  self._rna_scale  = self._robust_scale(rna_tr,  self.rna_min_scale)
        self._meth_center, self._meth_scale = self._robust_scale(meth_tr, self.meth_min_scale)
        self._cnv_center,  self._cnv_scale  = self._robust_scale(cnv_tr,  self.cnv_min_scale)

        # ------------------------------------------------------------------
        # Exposure standardisation — fit on training patients only.
        # Raw pack-years (0-160) must be standardised before entering the
        # ExposureEdgeBias and ExposureConditioner MLPs.  Without this,
        # large values (e.g. 157 pack-years) saturate the tanh in the
        # conditioner and produce near-identical outputs for all heavy
        # smokers, losing dose-response resolution above ~60 pack-years.
        # We use one row per patient (tumor row) to avoid counting each
        # patient twice.  std floor = 1.0 prevents division by zero if
        # all training patients have the same exposure value.
        # ------------------------------------------------------------------
        train_patient_meta = train_meta.drop_duplicates(
            subset=[self.patient_id_col]
        )
        exp_vals = train_patient_meta[self.exposure_cols].astype(
            np.float32
        ).to_numpy()                                    # [n_train, exp_dim]
        self._exp_mean = exp_vals.mean(axis=0)          # [exp_dim]
        self._exp_std  = np.maximum(
            exp_vals.std(axis=0), 1.0
        )                                               # [exp_dim], floor=1.0

        if self.verbose:
            for j, col in enumerate(self.exposure_cols):
                print(f"  Exposure [{col}]: "
                      f"train mean={self._exp_mean[j]:.2f}  "
                      f"train std={self._exp_std[j]:.2f}  "
                      f"→ standardised before model entry")

        self.scalers_fitted = True

        if self.verbose:
            print(f"[fit_scalers] {len(train_samples)} samples "
                  f"from {len(train_pids)} patients")
            floors = {"RNA": self.rna_min_scale,
                      "Meth": self.meth_min_scale,
                      "CNV":  self.cnv_min_scale}
            for name, raw, sc in [
                ("RNA",  rna_tr,  self._rna_scale),
                ("Meth", meth_tr, self._meth_scale),
                ("CNV",  cnv_tr,  self._cnv_scale),
            ]:
                pct = (sc == floors[name]).mean() * 100
                print(f"  {name:4s} | raw median={np.median(raw):.3g}  "
                      f"IQR median={np.median(sc):.3g}  "
                      f"IQR max={sc.max():.3g}  "
                      f"floor-clamped={pct:.0f}%")
                if name == "RNA" and np.median(sc) > 10:
                    print(f"  *** WARNING: RNA IQR={np.median(sc):.1f}. "
                          f"Apply log1p in Step_01. ***")
                # Thresholds are modality-specific:
                #   RNA  : >80% clamped = real problem (most genes should vary)
                #   Meth : >80% clamped at 0.005 floor = expected (promoters stable)
                #   CNV  : >99% clamped = problem; 90-99% is normal (most genes diploid)
                warn_threshold = {"RNA": 80, "Meth": 80, "CNV": 99}[name]
                if pct > warn_threshold:
                    print(f"  *** WARNING: {name} {pct:.0f}% floor-clamped "
                          f"(threshold {warn_threshold}%). "
                          f"Verify Step_01 output. ***")

    def _scale(self, x: np.ndarray,
               center: np.ndarray, scale: np.ndarray) -> np.ndarray:
        return np.tanh((x - center) / (scale * self.tanh_bound)).astype(np.float32)

    # ------------------------------------------------------------------
    # Graph construction
    # ------------------------------------------------------------------

    def _build_graph(self, sample_id: str) -> Data:
        if not self.scalers_fitted:
            raise RuntimeError(
                "Call dataset.fit_scalers(train_patient_ids) before iterating."
            )
        if sample_id not in self.rna.index:
            raise KeyError(f"Sample '{sample_id}' not in feature matrices.")

        rna_raw  = self.rna.loc[sample_id].to_numpy(dtype=np.float32)
        meth_raw = self.meth.loc[sample_id].to_numpy(dtype=np.float32)
        cnv_raw  = np.nan_to_num(
            self.cnv.loc[sample_id].to_numpy(dtype=np.float32), nan=0.0
        )
        cpg_raw  = self.cpg.loc[sample_id].to_numpy(dtype=np.float32)

        rna_s  = self._scale(rna_raw,  self._rna_center,  self._rna_scale)
        meth_s = self._scale(meth_raw, self._meth_center, self._meth_scale)
        cnv_s  = self._scale(cnv_raw,  self._cnv_center,  self._cnv_scale)
        # CpG beta is already in [0, 1]; NaN probes imputed as 0.5 (uncertain)
        cpg_s  = np.nan_to_num(cpg_raw, nan=0.5).astype(np.float32)

        # Gene nodes: [RNA, meth, CNV]  shape [n_genes, 3]
        gene_x = np.stack([rna_s, meth_s, cnv_s], axis=1)

        # CpG nodes: [beta, 0, 0]  shape [n_cpgs, 3]
        # Step_03 NodeTypeProjection uses only column 0 for CpG nodes
        cpg_x       = np.zeros((self.n_cpgs, 3), dtype=np.float32)
        cpg_x[:, 0] = cpg_s

        x = np.vstack([gene_x, cpg_x])   # [N, 3]

        return Data(
            x          = torch.tensor(x, dtype=torch.float32),
            node_type  = self.node_type,
            edge_index = self.edge_index,
            edge_attr  = self.edge_attr,
            edge_type  = self.edge_type,   # [E] 0=PPI, 1=CpG-Gene
        )

    # ------------------------------------------------------------------
    # Dataset interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.patient_ids)

    def __getitem__(self, idx: int):
        pid    = str(self.patient_ids[idx])
        rows   = self.meta[self.meta[self.patient_id_col] == pid]
        normal = rows[rows["_tissue"] == "normal"].iloc[0]
        tumor  = rows[rows["_tissue"] == "tumor"].iloc[0]

        normal_graph = self._build_graph(str(normal[self.sample_id_col]))
        tumor_graph  = self._build_graph(str(tumor[self.sample_id_col]))

        # Standardise exposure using training-set mean and std.
        # This prevents tanh saturation in ExposureEdgeBias and
        # ExposureConditioner MLPs for high pack-years patients.
        raw_exp = tumor[self.exposure_cols].astype(np.float32).to_numpy()
        exposure = torch.tensor(
            (raw_exp - self._exp_mean) / self._exp_std,
            dtype=torch.float32,
        )  # [exp_dim]  — zero-centred, unit-scale on training distribution

        return normal_graph, tumor_graph, exposure, pid

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def get_node_labels(self) -> Dict[int, str]:
        """Returns dict: node index -> human-readable label (gene or CpG)."""
        gene_map = pd.read_csv(os.path.join(self.data_dir, "map_gene_nodes.csv"))
        cpg_map  = pd.read_csv(os.path.join(self.data_dir, "map_cpg_nodes.csv"))
        labels   = {}
        for _, row in gene_map.iterrows():
            labels[int(row["node_idx"])] = f"GENE:{row['entrez_id']}"
        for _, row in cpg_map.iterrows():
            labels[int(row["node_idx"])] = f"CpG:{row['cpg_id']}"
        return labels

    def _print_summary(self) -> None:
        print("MultiOmicsDataset loaded:")
        print(f"  patients : {len(self.patient_ids)}")
        print(f"  genes    : {self.n_genes}")
        print(f"  CpGs     : {self.n_cpgs}")
        print(f"  nodes    : {self.n_nodes}")
        print(f"  exposure : {self.exposure_cols}")
