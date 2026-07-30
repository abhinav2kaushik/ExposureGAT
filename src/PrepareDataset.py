#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PrepareDataset.py  — degree-aware exposure gate

"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional
from torch import Tensor
from torch_geometric.nn import GATv2Conv
from torch_geometric.utils import softmax, degree


# =============================================================================
# GATv2Conv with direct logit bias injection  — UNCHANGED
# =============================================================================

class ExposureBiasGATv2Conv(GATv2Conv):
    """
    GATv2Conv subclass that injects a precomputed per-edge per-head bias
    directly into the scalar attention logit before softmax.
    UNCHANGED from previous version.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._logit_bias: Optional[Tensor] = None
        self.inspector.inspect_signature(self.edge_update)

    def set_logit_bias(self, bias: Optional[Tensor]) -> None:
        self._logit_bias = bias

    def edge_update(
        self,
        x_j:       Tensor,
        x_i:       Tensor,
        edge_attr: Optional[Tensor],
        index:     Tensor,
        ptr:       Optional[Tensor],
        dim_size:  Optional[int],
    ) -> Tensor:
        x = x_i + x_j

        if edge_attr is not None:
            if edge_attr.dim() == 1:
                edge_attr = edge_attr.view(-1, 1)
            assert self.lin_edge is not None
            edge_attr = self.lin_edge(edge_attr)
            edge_attr = edge_attr.view(-1, self.heads, self.out_channels)
            x = x + edge_attr

        x     = F.leaky_relu(x, self.negative_slope)
        alpha = (x * self.att).sum(dim=-1)              # [E, heads]

        if self._logit_bias is not None:
            alpha = alpha + self._logit_bias            # [E, heads]

        alpha = softmax(alpha, index, ptr, dim_size)
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)
        return alpha


# =============================================================================
# Learnable soft-threshold degree transform  — UNCHANGED
# =============================================================================

class DegreeTransform(nn.Module):
    """
    Learnable sigmoid that maps raw node degree to (0, 1).

    degree << threshold  ->  ~0   peripheral nodes suppressed
    degree ~  threshold  ->  0.5  transition zone
    degree >> threshold  ->  ~1   hub nodes, saturated (all equal)

    Initialised from biological prior (threshold=20, steepness=0.1).
    Clamped during training: steepness in [0.01, 0.5], threshold in [1, 200].
    """

    def __init__(self, init_threshold: float = 20.0,
                       init_steepness: float = 0.1):
        super().__init__()
        self.threshold = nn.Parameter(torch.tensor(init_threshold))
        self.steepness = nn.Parameter(torch.tensor(init_steepness))

    def forward(self, deg: Tensor) -> Tensor:
        s = torch.clamp(self.steepness, min=0.01, max=0.5)
        t = torch.clamp(self.threshold, min=1.0,  max=200.0)
        return torch.sigmoid(s * (deg - t))              # (0, 1)


# =============================================================================
# ExposureEdgeBias — FIX 3 + FIX 4 applied
# =============================================================================

class ExposureEdgeBias(nn.Module):
    """
    Maps (exposure, edge_type, src_degree, tgt_degree) to a per-edge
    per-head attention logit bias.

    Architecture: two separate pathways combined additively.

        Pathway 1 — Main gate (smoking DRIVES this):
            [smoking | edge_type_one_hot] -> MLP -> main_logit

        Pathway 2 — Degree nudge (degree TWEAKS, does NOT drive):
            [src_deg_transformed | tgt_deg_transformed]
            -> Linear (no hidden, no ReLU)
            -> tanh()                      <- FIX 4: hard bounds nudge
            -> * nudge_scale               <- true cap at [-nudge_scale, +nudge_scale]

        Final:  (FIX 3: was sigmoid, now 0.5*tanh for signed output)
            gate = 0.5 * tanh(main_logit + nudge_scale * tanh(nudge))

    FIX 3 — signed bias:
        Previous: sigmoid(logit) -> always (0,1) -> gate only adds positive bias.
        Fixed:    0.5 * tanh(logit) -> (-0.5, +0.5) -> smoking can activate OR
                  suppress attention. Biologically: smoking may silence some
                  regulatory edges while activating others.

    FIX 4 — bounded nudge:
        Previous: nudge_scale * nudge, where nudge is unbounded (trainable weights).
        Fixed:    nudge_scale * tanh(nudge) -> nudge strictly in [-nudge_scale, +nudge_scale].
                  Degree cannot escape its intended role as a small correction
                  regardless of how large degree_nudge weights grow.
    """

    N_EDGE_TYPES = 2   # 0 = PPI, 1 = CpG-Gene

    def __init__(self, exp_dim: int, heads: int,
                 hidden: int = 32, nudge_scale: float = 0.1):
        super().__init__()

        self.nudge_scale   = nudge_scale
        self.deg_transform = DegreeTransform(
            init_threshold = 20.0,
            init_steepness = 0.1,
        )

        # Pathway 1: main gate — smoking drives attention modulation
        self.main_gate = nn.Sequential(
            nn.Linear(exp_dim + self.N_EDGE_TYPES, hidden),
            nn.ReLU(),
            nn.Linear(hidden, heads),   # raw logit, NOT activated yet
        )

        # Pathway 2: degree nudge — low capacity, near-zero init
        self.degree_nudge = nn.Linear(2, heads, bias=False)
        nn.init.normal_(self.degree_nudge.weight, mean=0.0, std=0.01)

    def forward(
        self,
        exposure:   Tensor,   # [B, exp_dim]
        edge_type:  Tensor,   # [E]
        edge_batch: Tensor,   # [E]
        src_deg:    Tensor,   # [E]  float, already local-indexed
        tgt_deg:    Tensor,   # [E]  float, already local-indexed
    ) -> Tensor:
        """Returns [E, heads] logit bias. Output range: (-0.5, +0.5) per head."""

        # --- Pathway 1: main logit ---
        etype_oh   = F.one_hot(edge_type.long(),
                               num_classes=self.N_EDGE_TYPES).float()  # [E, 2]
        main_input = torch.cat([exposure[edge_batch], etype_oh], dim=1) # [E, exp_dim+2]
        main_logit = self.main_gate(main_input)                         # [E, heads]

        # --- Pathway 2: degree nudge ---
        # DETACH: structural context, not a gradient pathway
        src_d = self.deg_transform(
            src_deg.float().detach()
        ).unsqueeze(-1)                                                 # [E, 1]
        tgt_d = self.deg_transform(
            tgt_deg.float().detach()
        ).unsqueeze(-1)                                                 # [E, 1]

        degree_input = torch.cat([src_d, tgt_d], dim=-1)               # [E, 2]

        # FIX 4: tanh hard-bounds the nudge before scaling
        nudge = torch.tanh(self.degree_nudge(degree_input))             # [E, heads] in (-1,+1)

        # --- Combine ---
        final_logit = main_logit + self.nudge_scale * nudge             # [E, heads]

        # FIX 3: 0.5*tanh gives signed output in (-0.5, +0.5)
        # Smoking can now SUPPRESS edges (negative bias) or ACTIVATE them
        return 0.5 * torch.tanh(final_logit)                            # [E, heads]


# =============================================================================
# NodeTypeProjection  — UNCHANGED
# =============================================================================

class NodeTypeProjection(nn.Module):
    def __init__(self, cpg_in: int, gene_in: int, hidden_dim: int):
        super().__init__()
        self.cpg_proj  = nn.Linear(cpg_in,  hidden_dim)
        self.gene_proj = nn.Linear(gene_in, hidden_dim)

    def forward(self, x: Tensor, node_type: Tensor) -> Tensor:
        out = torch.zeros(x.size(0), self.cpg_proj.out_features,
                          device=x.device, dtype=x.dtype)
        cpg_mask  = node_type == 0
        gene_mask = node_type == 1
        if cpg_mask.any():
            out[cpg_mask]  = self.cpg_proj(x[cpg_mask, :1])
        if gene_mask.any():
            out[gene_mask] = self.gene_proj(x[gene_mask, :3])
        return out


# =============================================================================
# ExposureConditioner  — UNCHANGED
# =============================================================================

class ExposureConditioner(nn.Module):
    def __init__(self, exp_dim: int, out_dim: int,
                 hidden: int = 64, tanh_scale: float = 0.5):
        super().__init__()
        self.tanh_scale = tanh_scale
        self.mlp = nn.Sequential(
            nn.Linear(exp_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, out_dim),
        )
        nn.init.zeros_(self.mlp[-1].weight)
        nn.init.zeros_(self.mlp[-1].bias)

    def forward(self, h: Tensor, exposure: Tensor,
                node_batch: Tensor) -> Tensor:
        scale = self.mlp(exposure)
        return h * (1.0 + self.tanh_scale * torch.tanh(scale[node_batch]))


# =============================================================================
# ExposureGAT — FIX 1 + FIX 5 applied
# =============================================================================

class ExposureGAT(nn.Module):
    """
    Exposure-conditioned GAT with degree-aware gate, all bugs fixed.

    FIX 1 — Batch indexing:
        PyG offsets node indices in batched graphs.
        edge_index[0] for batch of size B contains indices 0 .. B*N-1.
        node_degree is only length N.
        Fix: local_idx = batched_idx % self.n_nodes before indexing.
        Safe because all patients share identical topology and node order.

    FIX 5 — Directional degree:
        Two separate degree buffers precomputed at construction:
          node_out_degree[i] = number of outgoing edges from node i
          node_in_degree[i]  = number of incoming edges to   node i
        source nodes use out-degree (how many signals they broadcast)
        target nodes use in-degree  (how many signals they receive;
          directly relevant to softmax competition over incoming edges)
        For symmetric/bidirectional graphs in/out degrees are equal,
        so results are identical — but the semantics are now correct.
    """

    def __init__(
        self,
        cpg_in:      int   = 1,
        gene_in:     int   = 3,
        exp_dim:     int   = 1,
        hidden_dim:  int   = 64,
        out_dim:     int   = 32,
        heads:       int   = 4,
        dropout:     float = 0.1,
        tanh_scale:  float = 0.5,
        bias_hidden: int   = 32,
        nudge_scale: float = 0.1,
        n_nodes:               int            = 0,
        edge_index_for_degree: Optional[Tensor] = None,
    ):
        super().__init__()
        self.out_dim = out_dim
        self.heads   = heads

        # ------------------------------------------------------------------
        # Guard: n_nodes and edge_index_for_degree are mandatory.
        # Silent fallback to zeros is worse than a clear error — it trains
        # a model where the degree gate is silently disabled, with no
        # indication that anything is wrong.
        # ------------------------------------------------------------------
        if n_nodes <= 0:
            raise ValueError(
                "ExposureGAT requires n_nodes > 0 for degree-aware batching. "
                "Pass n_nodes=dataset.n_nodes from Step_04/05."
            )
        if edge_index_for_degree is None:
            raise ValueError(
                "ExposureGAT requires edge_index_for_degree for degree "
                "precomputation. Pass edge_index_for_degree=dataset.edge_index "
                "from Step_04/05."
            )

        self.n_nodes = n_nodes           # FIX 1: stored for modulo indexing

        # Input projection
        self.input_proj = NodeTypeProjection(cpg_in, gene_in, hidden_dim)

        # GAT layers
        self.conv1 = ExposureBiasGATv2Conv(
            in_channels    = hidden_dim,
            out_channels   = hidden_dim,
            heads          = heads,
            edge_dim       = 1,
            dropout        = dropout,
            add_self_loops = False,
        )
        self.conv2 = ExposureBiasGATv2Conv(
            in_channels    = hidden_dim * heads,
            out_channels   = out_dim,
            heads          = 1,
            edge_dim       = 1,
            dropout        = dropout,
            add_self_loops = False,
        )

        # Edge bias modules
        self.edge_bias1 = ExposureEdgeBias(
            exp_dim=exp_dim, heads=heads,
            hidden=bias_hidden, nudge_scale=nudge_scale,
        )
        self.edge_bias2 = ExposureEdgeBias(
            exp_dim=exp_dim, heads=1,
            hidden=bias_hidden, nudge_scale=nudge_scale,
        )

        # Node conditioner
        self.conditioner = ExposureConditioner(
            exp_dim, out_dim, hidden=64, tanh_scale=tanh_scale
        )

        # Layer norms
        self.norm1 = nn.LayerNorm(hidden_dim * heads)
        self.norm2 = nn.LayerNorm(out_dim)

        # ------------------------------------------------------------------
        # FIX 5: precompute BOTH out-degree and in-degree separately.
        #
        # out_degree[i] = |{j : (i,j) in E}|  — signals node i broadcasts
        # in_degree[i]  = |{j : (j,i) in E}|  — signals node i receives
        #
        # Stored as non-trainable buffers: move to GPU with model.to(device),
        # saved in state dict, never updated by optimizer.
        # ------------------------------------------------------------------
        # Guards passed above — edge_index_for_degree is not None and n_nodes > 0
        out_deg = degree(
            edge_index_for_degree[0],
            num_nodes = n_nodes,
            dtype     = torch.float32,
        )
        in_deg = degree(
            edge_index_for_degree[1],
            num_nodes = n_nodes,
            dtype     = torch.float32,
        )

        self.register_buffer("node_out_degree", out_deg)   # [N]
        self.register_buffer("node_in_degree",  in_deg)    # [N]

    def _encode(self, data, exposure: Tensor):
        x          = torch.nan_to_num(data.x, nan=0.0)
        node_type  = data.node_type
        edge_index = data.edge_index
        edge_attr  = data.edge_attr.view(-1, 1)
        edge_type  = data.edge_type
        batch      = data.batch
        edge_batch = batch[edge_index[0]]

        # ------------------------------------------------------------------
        # FIX 1 — convert batched node indices to local graph indices
        # before looking up node_degree.
        #
        # PyG batching adds k*N to all node indices of graph k.
        # node_out_degree and node_in_degree are length N (single graph).
        # Modulo maps batched indices back to [0, N-1].
        # Safe because all patients share identical node ordering.
        # ------------------------------------------------------------------
        local_src = edge_index[0] % self.n_nodes    # [E] in [0, N-1]
        local_tgt = edge_index[1] % self.n_nodes    # [E] in [0, N-1]

        # FIX 5 — source uses out-degree, target uses in-degree
        src_deg = self.node_out_degree[local_src]   # [E]
        tgt_deg = self.node_in_degree[local_tgt]    # [E]

        # 1. Input projection
        h = F.elu(self.input_proj(x, node_type))    # [N_batch, hidden_dim]

        # 2. Conv1
        self.conv1.set_logit_bias(
            self.edge_bias1(exposure, edge_type, edge_batch, src_deg, tgt_deg)
        )
        h1, (ei, a1h) = self.conv1(
            h, edge_index, edge_attr, return_attention_weights=True
        )
        alpha1 = a1h.mean(dim=1)
        h1     = self.norm1(F.elu(h1))

        # 3. Conv2
        self.conv2.set_logit_bias(
            self.edge_bias2(exposure, edge_type, edge_batch, src_deg, tgt_deg)
        )
        h2, (_, a2h) = self.conv2(
            h1, edge_index, edge_attr, return_attention_weights=True
        )
        alpha2 = a2h.squeeze(1)
        h2     = self.norm2(F.elu(h2))

        alpha = alpha1 + alpha2
        h_out = self.conditioner(h2, exposure, batch)

        return h_out, alpha, ei

    def forward(self, normal, tumor, exposure: Tensor) -> dict:
        """UNCHANGED. Step_04 and Step_05 require only constructor changes."""
        baseline = torch.zeros_like(exposure)

        h_n, alpha_n, edge_index = self._encode(normal, baseline)
        h_t, alpha_t, _          = self._encode(tumor,  exposure)

        return {
            "h_normal":    h_n,
            "h_tumor":     h_t,
            "alpha_n":     alpha_n,
            "alpha_t":     alpha_t,
            "delta_h":     h_t - h_n,
            "delta_alpha": alpha_t - alpha_n,
            "edge_index":  edge_index,
        }

    def report_degree_thresholds(self) -> None:
        """
        Print the learned degree transform parameters for both GAT layers.

        Call after training to inspect what the model learned:
            layer-1 threshold: degree cutoff in low-level structural features
            layer-2 threshold: degree cutoff in high-level aggregated features

        These may differ because layer 1 and layer 2 capture different levels
        of abstraction. Do not report a single learned hub threshold in results
        without inspecting both values.

        Example interpretation:
            layer-1 threshold = 8   -> even low-degree nodes matter early
            layer-2 threshold = 35  -> only true hubs matter at the output level
        """
        for layer_name, bias_module in [("layer-1 (conv1)", self.edge_bias1),
                                         ("layer-2 (conv2)", self.edge_bias2)]:
            t = bias_module.deg_transform.threshold.item()
            s = bias_module.deg_transform.steepness.item()
            s_clamped = max(0.01, min(0.5,  s))
            t_clamped = max(1.0,  min(200.0, t))
            print(f"  DegreeTransform {layer_name}:")
            print(f"    threshold (raw={t:.2f}, clamped={t_clamped:.2f})"
                  f"  — degree above this is treated as a hub")
            print(f"    steepness (raw={s:.4f}, clamped={s_clamped:.4f})"
                  f"  — sharpness of the hub transition")
            # Quick sanity: show transform values at key degree points
            import torch as _torch
            probe = _torch.tensor([1., 5., 10., 20., 50., 100., 200.])
            vals  = _torch.sigmoid(
                s_clamped * (probe - t_clamped)
            ).tolist()
            print(f"    transform values at degree "
                  f"[1,5,10,20,50,100,200]: "
                  f"{[round(v,3) for v in vals]}")
        print(f"  nudge_scale: {self.edge_bias1.nudge_scale}  "
              f"(max logit shift from degree = ±{self.edge_bias1.nudge_scale})")
