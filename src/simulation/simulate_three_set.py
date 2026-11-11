#!/usr/bin/env python3
"""Exposure-associated redistribution benchmark. Requires numpy, pandas, scipy.
Simulation used in the manuscript (Figure 2). Set A: paired changes redistributed according to
exposure with the total planted change held constant; Set C: the same shifts assigned independently
of exposure; Set B: background. Regular graph (equal degree for all nodes). See README.md.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata


def regular_graph(n, degree, rng):
    if degree < 2 or degree >= n or degree % 2:
        raise ValueError('degree must be even and satisfy 2 <= degree < n')
    edges = {tuple(sorted((i, (i + j) % n))) for i in range(n) for j in range(1, degree // 2 + 1)}
    edges = list(sorted(edges))
    present = set(edges)
    for _ in range(20 * len(edges)):
        i, j = rng.choice(len(edges), 2, replace=False)
        a, b = edges[i]; c, d = edges[j]
        if rng.random() < .5:
            c, d = d, c
        if len({a, b, c, d}) != 4:
            continue
        x, y = tuple(sorted((a, d))), tuple(sorted((c, b)))
        if x in present or y in present:
            continue
        present.remove(edges[i]); present.remove(edges[j])
        edges[i], edges[j] = x, y
        present.update((x, y))
    return np.array(sorted(present), dtype=int)


def plant(n_items, n_a, n_c, exposure, rng, null=False):
    if min(n_a, n_c) < 2 or n_a % 2 or n_c % 2 or n_a + n_c >= n_items:
        raise ValueError('A and C counts must be positive even numbers with A+C < item count')
    order = rng.permutation(n_items)
    A, C = order[:n_a], order[n_a:n_a+n_c]
    labels = np.full(n_items, 'B'); labels[A] = 'A'; labels[C] = 'C'
    thresholds = rng.uniform(.15, .70, n_a // 2)
    dose = exposure / (exposure + 40.)
    swap = dose[:, None] > thresholds[None, :]
    if null:
        # Independent latent variation retains A heterogeneity without exposure dependence.
        swap = rng.random((len(exposure), n_a // 2)) < .5
    delta = np.zeros((len(exposure), n_items))
    low, high = .4, 1.6
    delta[:, A[::2]] = np.where(swap, high, low)
    delta[:, A[1::2]] = np.where(swap, low, high)
    delta[:, C[::2]] = low; delta[:, C[1::2]] = high
    return delta, labels, A.reshape(-1, 2), thresholds


def auc(y, score):
    y = np.asarray(y, bool)
    n1, n0 = y.sum(), (~y).sum()
    return float((rankdata(score)[y].sum() - n1*(n1+1)/2)/(n1*n0))


def corr(x, y):
    if np.ptp(x) < 1e-10 or np.ptp(y) < 1e-10:
        return 0.
    return float(spearmanr(x, y).statistic)


def diagnostics(delta, observed, labels, exposure):
    scores = np.array([abs(corr(observed[:, j], exposure)) for j in range(observed.shape[1])])
    ac = labels != 'B'
    return {'latent_L1_range': float(np.ptp(np.abs(delta).sum(1))),
            'latent_L2_squared_range': float(np.ptp((delta**2).sum(1))),
            'observed_L1_exposure_rho': corr(np.abs(observed).sum(1), exposure),
            'observed_L2_squared_exposure_rho': corr((observed**2).sum(1), exposure),
            'A_vs_C_association_AUROC': auc(labels[ac] == 'A', scores[ac]),
            'A_vs_B_association_AUROC': auc(labels[labels != 'C'] == 'A', scores[labels != 'C']),
            'C_vs_B_mean_change_AUROC': auc(labels[labels != 'A'] == 'C', np.abs(observed.mean(0))[labels != 'A'])}


def generate(out_dir, n_patients=250, n_genes=400, n_a=80, n_c=80,
             degree=8, effect=1., sigma=.7, seed=1, null=False):
    if n_patients < 10 or effect < 0 or sigma <= 0:
        raise ValueError('Require n_patients >= 10, effect >= 0, sigma > 0')
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    exposure = np.minimum(rng.exponential(39., n_patients), 160.)
    exposure[rng.random(n_patients) < .12] = 0
    und = regular_graph(n_genes, degree, rng)
    latent, labels, pairs, thresholds = plant(n_genes, n_a, n_c, exposure, rng, null)
    latent *= effect
    patients = [f'P{i+1:04d}' for i in range(n_patients)]
    samples = [f'{p}-{t}' for p in patients for t in ('NT', 'TP')]
    genes = [str(100000+i) for i in range(n_genes)]
    def write_features(name, normal, tumor):
        data = np.empty((2*n_patients, n_genes))
        data[::2], data[1::2] = normal, tumor
        pd.DataFrame(data, index=samples, columns=genes).to_csv(out/name)
    results = {}
    # Exchangeable measurement noise is essential to the norm-invariance guarantee.
    for modality, scale, noise in [('rna', sigma, sigma/np.sqrt(2)), ('meth', -.06, .02), ('cnv', .10, .04)]:
        shape = (n_patients, n_genes)
        if modality == 'meth':
            if effect > 3:
                raise ValueError('effect > 3 would violate bounded methylation construction')
            baseline = .5 + rng.uniform(-.05, .05, shape)
            normal = baseline + rng.uniform(-noise, noise, shape)
            tumor = baseline + rng.uniform(-noise, noise, shape) + scale*latent
            assert normal.min() >= 0 and tumor.min() >= 0 and normal.max() <= 1 and tumor.max() <= 1
        else:
            baseline = (4.2 if modality == 'rna' else 0.) + rng.normal(0, .3, shape)
            normal = baseline + rng.normal(0, noise, shape)
            tumor = baseline + rng.normal(0, noise, shape) + scale*latent
        write_features(f'X_{modality}_gene_features.csv', normal, tumor)
        results[modality] = diagnostics(scale*latent, tumor-normal, labels, exposure)
    pd.DataFrame({'sample_id': samples, 'patient_id': np.repeat(patients, 2),
                  'tissue_type': ['NT','TP']*n_patients,
                  'smoking_intensity': np.repeat(exposure, 2)}).to_csv(out/'metadata_samples.csv', index=False)
    pd.DataFrame(np.vstack([und, und[:, ::-1]]), columns=['from_idx','to_idx']).to_csv(out/'edges_ppi.csv', index=False)
    pd.DataFrame({'entrez_id': genes, 'node_idx': np.arange(n_genes)}).to_csv(out/'map_gene_nodes.csv', index=False)
    # Gene-level methylation is included; no additional CpG nodes in this benchmark.
    pd.DataFrame(index=samples).to_csv(out/'X_cpg_features.csv')
    pd.DataFrame(columns=['cpg_id','node_idx']).to_csv(out/'map_cpg_nodes.csv', index=False)
    pd.DataFrame(columns=['from_idx','to_idx']).to_csv(out/'edges_cpg_to_gene.csv', index=False)
    pd.DataFrame({'node_idx': np.arange(n_genes), 'node_label': ['GENE:'+g for g in genes],
                  'node_type':'gene', 'set_membership': labels, 'degree': degree,
                  'exposure_truth': (labels == 'A') & (not null)}).to_csv(out/'sim_ground_truth_nodes.csv', index=False)
    pd.DataFrame({'pair_node_1':pairs[:,0], 'pair_node_2':pairs[:,1], 'dose_threshold':thresholds}).to_csv(out/'sim_A_pairs.csv',index=False)
    np.savez_compressed(out/'sim_latent_node_changes.npz', delta=latent)
    # Independent edge-level task, not inferred merely from touching an A node.
    edge_delta, edge_labels, _, _ = plant(len(und), n_a, n_c, exposure, rng, null)
    edge_delta *= .1*effect
    edge_n = 1 + rng.normal(0,.02,edge_delta.shape)
    edge_t = 1 + rng.normal(0,.02,edge_delta.shape) + edge_delta
    pd.DataFrame({'edge_id': np.arange(len(und)), 'from_idx':und[:,0], 'to_idx':und[:,1],
                  'set_membership':edge_labels, 'exposure_truth':(edge_labels=='A') & (not null)}).to_csv(out/'sim_ground_truth_edges.csv',index=False)
    pd.DataFrame(edge_n,index=patients).to_csv(out/'edge_weights_normal.csv')
    pd.DataFrame(edge_t,index=patients).to_csv(out/'edge_weights_tumor.csv')
    np.savez_compressed(out/'sim_latent_edge_changes.npz',delta=edge_delta)
    results['edges'] = diagnostics(edge_delta,edge_t-edge_n,edge_labels,exposure)
    degree_observed = np.bincount(und.ravel(),minlength=n_genes)
    assert np.all(degree_observed == degree)
    for item in results.values():
        assert item['latent_L1_range'] < 1e-8 and item['latent_L2_squared_range'] < 1e-8
    config = dict(n_patients=n_patients,n_genes=n_genes,n_a=n_a,n_c=n_c,degree=degree,effect=effect,sigma=sigma,seed=seed,null=null)
    (out/'simulation_config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'design_checks.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out_dir', required=True)
    for name, default in [('n_patients',250),('n_genes',400),('n_a',80),('n_c',80),('degree',8),('seed',1)]:
        p.add_argument('--'+name,type=int,default=default)
    p.add_argument('--effect',type=float,default=1.)
    p.add_argument('--sigma',type=float,default=.7)
    p.add_argument('--null',action='store_true')
    print(json.dumps(generate(**vars(p.parse_args())),indent=2))
if __name__ == '__main__':
    main()
