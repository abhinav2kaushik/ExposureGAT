"""Simulate a small paired normal/tumor multi-omics dataset in ExposureGAT input format.

Used for the demo run and for testing. The design plants three node sets:
  A  exposure-linked: tumor-minus-normal change scales with patient exposure
  C  fixed transition: the same tumor shift in every patient, unrelated to exposure
  B  background: noise only
A recovers well when the model ranks exposure-linked change above generic
tumor-normal change (A vs C is the informative comparison).

Output files match what src/data_loader.py expects, plus sim_ground_truth_nodes.csv.

Example:
  python src/simulation/simulate_data.py --out_dir tmp/demo_data --seed 1
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def scale_free_edges(n, m, rng):
    """Undirected preferential-attachment graph on n nodes (each new node adds m edges)."""
    edges, targets, repeated = set(), list(range(m)), []
    for new in range(m, n):
        for t in set(targets):
            edges.add((min(new, t), max(new, t)))
        repeated.extend(targets + [new] * m)
        targets = list(rng.choice(repeated, size=m, replace=False))
    return sorted(edges)


def main(a):
    rng = np.random.default_rng(a.seed)
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output folder must be new or empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    P, G, C = a.patients, a.genes, a.cpgs

    # --- exposure: pack-years, with confirmed never-smokers at zero -------------
    never = rng.random(P) < a.never_fraction
    pack_years = np.where(never, 0.0, np.minimum(rng.gamma(2.0, 20.0, P), 150.0)).round(1)
    status = np.where(never, 'Lifelong Non-Smoker',
                      rng.choice(['Current Smoker', 'Current Reformed Smoker for < or = 15 yrs',
                                  'Current Reformed Smoker for > 15 yrs'], P))
    z = (pack_years - pack_years.mean()) / pack_years.std()
    age = rng.integers(45, 85, P)
    sex = rng.choice(['male', 'female'], P)

    # --- node sets ----------------------------------------------------------------
    order = rng.permutation(G)
    set_a, set_c = order[:a.n_exposure_genes], order[a.n_exposure_genes:a.n_exposure_genes + a.n_transition_genes]
    membership = np.array(['B'] * G, dtype=object)
    membership[set_a], membership[set_c] = 'A', 'C'

    # --- graph: gene-gene (scale-free) and CpG -> nearest gene ---------------------
    ppi = scale_free_edges(G, a.edges_per_gene, rng)
    cpg_gene = rng.integers(0, G, C)
    linked_a = np.isin(cpg_gene, set_a)
    cpg_membership = np.where(linked_a, 'A', 'B')

    # --- molecular profiles ---------------------------------------------------------
    gene_mean = rng.normal(6.0, 1.5, G)
    meth_mean = rng.uniform(0.1, 0.9, G)
    cpg_mean = rng.uniform(0.1, 0.9, C)
    rows, rna, meth, cnv, cpg = [], [], [], [], []
    for p in range(P):
        patient_rna = gene_mean + rng.normal(0, 0.4, G)
        patient_meth = np.clip(meth_mean + rng.normal(0, 0.03, G), 0.01, 0.99)
        patient_cpg = np.clip(cpg_mean + rng.normal(0, 0.03, C), 0.01, 0.99)
        for tissue in ('NT', 'TP'):
            r = patient_rna + rng.normal(0, 0.25, G)
            m = np.clip(patient_meth + rng.normal(0, 0.02, G), 0, 1)
            c = np.zeros(G)
            b = np.clip(patient_cpg + rng.normal(0, 0.02, C), 0, 1)
            if tissue == 'TP':
                r = r + rng.normal(0, 0.5, G)                       # generic tumor noise
                r[set_c] += a.transition_effect                    # same shift in all patients
                r[set_a] += a.exposure_effect * z[p] + rng.normal(0, 0.2, len(set_a))
                m[set_a] = np.clip(m[set_a] - 0.04 * z[p], 0, 1)   # exposure-linked methylation loss
                events = rng.random(G) < 0.08                      # tumor copy-number events
                c[events] = np.sign(rng.normal(size=events.sum())) * np.log2(rng.integers(1, 4, events.sum()) + 1)
                b = np.clip(b + rng.normal(0, 0.05, C), 0, 1)
                b[linked_a] = np.clip(b[linked_a] - 0.05 * z[p], 0, 1)
            sample = f'SIM{p:03d}-{tissue}'
            rows.append(dict(sample_id=sample, patient_id=f'SIM{p:03d}', tissue_type=tissue,
                             smoking_intensity=pack_years[p], smoking_status=status[p],
                             age=int(age[p]), gender=sex[p]))
            rna.append(r); meth.append(m); cnv.append(c); cpg.append(b)

    genes = [str(10000 + g) for g in range(G)]
    cpgs = [f'cgSIM{c:05d}' for c in range(C)]
    samples = [r['sample_id'] for r in rows]
    for name, values, cols in [('rna_gene', rna, genes), ('meth_gene', meth, genes),
                               ('cnv_gene', cnv, genes), ('cpg', cpg, cpgs)]:
        pd.DataFrame(np.array(values), index=samples, columns=cols).to_csv(out / f'X_{name}_features.csv')
    pd.DataFrame(rows).to_csv(out / 'metadata_samples.csv', index=False)

    both = [(u, v) for u, v in ppi] + [(v, u) for u, v in ppi]
    pd.DataFrame(both, columns=['from_idx', 'to_idx']).assign(weight=1.0).to_csv(out / 'edges_ppi.csv', index=False)
    pd.DataFrame(dict(from_idx=np.arange(C) + G, to_idx=cpg_gene, weight=1.0, cpg_id=cpgs,
                      entrez_id=[genes[g] for g in cpg_gene])).to_csv(out / 'edges_cpg_to_gene.csv', index=False)
    pd.DataFrame(dict(entrez_id=genes, node_idx=np.arange(G), node_type=1)).to_csv(out / 'map_gene_nodes.csv', index=False)
    pd.DataFrame(dict(cpg_id=cpgs, node_idx=np.arange(C) + G, node_type=0)).to_csv(out / 'map_cpg_nodes.csv', index=False)
    pd.DataFrame(dict(node_idx=np.arange(G + C), node_label=['GENE:' + g for g in genes] + ['CpG:' + c for c in cpgs],
                      set_membership=list(membership) + list(cpg_membership))).to_csv(out / 'sim_ground_truth_nodes.csv', index=False)
    print(f'Simulated {P} paired patients, {G} genes ({len(set_a)} exposure-linked, {len(set_c)} fixed-transition), '
          f'{C} CpGs, {2 * len(ppi) + 2 * C} directed edges -> {out}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out_dir', required=True)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--patients', type=int, default=60)
    p.add_argument('--genes', type=int, default=300)
    p.add_argument('--cpgs', type=int, default=30)
    p.add_argument('--edges_per_gene', type=int, default=3)
    p.add_argument('--n_exposure_genes', type=int, default=20)
    p.add_argument('--n_transition_genes', type=int, default=20)
    p.add_argument('--exposure_effect', type=float, default=0.8)
    p.add_argument('--transition_effect', type=float, default=1.0)
    p.add_argument('--never_fraction', type=float, default=0.15)
    main(p.parse_args())
