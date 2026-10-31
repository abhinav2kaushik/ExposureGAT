"""Simulation-only scaling adapter (selected with exposureGAT.py --scaling shared).

Pairing, exposure standardization, tanh, graph building, and training-only fitting
remain inherited unchanged. Do not automatically use this scaling for TCGA.
"""
import numpy as np
from data_loader import LungMultiOmicsDataset


class SimulationDataset(LungMultiOmicsDataset):
    @staticmethod
    def _robust_scale(X, min_scale):
        # Shared across genes, fitted separately per modality using training samples.
        center = float(np.median(X))
        q25, q75 = np.percentile(X, [25, 75])
        scale = max(float(q75-q25), min_scale)
        return (np.full(X.shape[1], center, dtype=np.float32),
                np.full(X.shape[1], scale, dtype=np.float32))
