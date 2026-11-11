#!/usr/bin/env python3
"""Measured direction of tumor-minus-normal change for the top genes by weighted node score.

RNA: paired difference (transition minus reference) in each patient; median across patients and
fraction of patients with a higher transition value. Promoter methylation: restricted to patients
whose reference and transition methylation profiles were both measured (profiles imputed by the
cohort normal-median are detected as identical across samples); genes without promoter-annotated
probes (value equal to the within-sample median in every sample) are reported as not measured.
The direction describes measured molecular change, not the sign of the model contribution.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main(a):
    data, out = Path(a.data), Path(a.out); out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(data / 'metadata_samples.csv')
    ref, tra = a.reference_label, a.transition_label
    pairs = meta.pivot_table(index='patient_id', columns='tissue_type', values='sample_id', aggfunc='first').dropna(subset=[ref, tra])
    rna = pd.read_csv(data / 'X_rna_gene_features.csv', index_col=0); rna.columns = rna.columns.astype(str)
    meth = pd.read_csv(data / 'X_meth_gene_features.csv', index_col=0); meth.columns = meth.columns.astype(str)
    cons = pd.read_csv(Path(a.weighted) / 'weighted_node_scores.csv', dtype={'entrez_id': str})
    top = cons[cons['rank'] <= int(np.ceil(a.top_fraction * len(cons)))]
    d_rna = pd.DataFrame(rna.loc[pairs[tra]].values - rna.loc[pairs[ref]].values, index=pairs.index, columns=rna.columns)
    imputed = meth.round(6).duplicated(keep=False)
    complete = [p for p in pairs.index if not imputed[pairs.loc[p, ref]] and not imputed[pairs.loc[p, tra]]]
    d_meth = pd.DataFrame(meth.loc[pairs.loc[complete, tra]].values - meth.loc[pairs.loc[complete, ref]].values,
                          index=complete, columns=meth.columns)
    no_probe = meth.sub(meth.median(axis=1), axis=0).abs().max() < 1e-9
    rows = []
    for r in top.itertuples():
        e = r.entrez_id
        x = d_rna[e] if e in d_rna else pd.Series(dtype=float)
        measured = e in d_meth and not no_probe.get(e, True) and len(complete) > 0
        m = d_meth[e] if measured else pd.Series(dtype=float)
        rows.append({'Rank (weighted node score)': r.rank, 'Entrez ID': e, 'Weighted node score': r.weighted_node_score,
                     'Exposure alignment (Spearman)': r.exposure_alignment,
                     'RNA: median transition−reference': x.median(), 'RNA: fraction of patients higher in transition': (x > 0).mean(),
                     'RNA direction': 'up' if x.median() > 0 else 'down',
                     'Methylation: complete pairs (n)': len(complete),
                     'Promoter methylation: median transition−reference': m.median() if measured else np.nan,
                     'Methylation: fraction higher in transition': (m > 0).mean() if measured else np.nan,
                     'Methylation direction': ('up' if m.median() > 0 else 'down') if measured else 'not measured'})
    t = pd.DataFrame(rows); t.to_csv(out / 'direction_of_change_top_genes.csv', index=False)
    print(f"{len(t)} genes; RNA {t['RNA direction'].value_counts().to_dict()}; methylation "
          f"{t['Methylation direction'].value_counts().to_dict()} ({len(complete)} complete pairs)")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data', required=True); p.add_argument('--weighted', required=True); p.add_argument('--out', required=True)
    p.add_argument('--top_fraction', type=float, default=.10)
    p.add_argument('--reference_label', default='NT'); p.add_argument('--transition_label', default='TP')
    main(p.parse_args())
