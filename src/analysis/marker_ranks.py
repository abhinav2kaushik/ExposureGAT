#!/usr/bin/env python3
"""Positions of prespecified marker genes and CpG probes (for example, smoking anchors) in the
ExposureGAT rankings: weighted node score rank and components (genes ranked among genes, CpGs among CpGs), cohort node-priority rank, and
whether each marker is present in the graph."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cons = pd.read_csv(Path(a.weighted) / 'weighted_node_scores.csv', dtype={'entrez_id': str}).set_index('entrez_id')
    cpg_file = Path(a.weighted) / 'weighted_node_scores_CpG.csv'
    cpg = pd.read_csv(cpg_file).set_index('node_label') if cpg_file.exists() else pd.DataFrame()
    nodes = pd.read_csv(Path(a.model) / 'cohort_summary.csv').set_index('node_label')
    pr = nodes.median_patient_rank.rank(method='min')
    rows = []
    for g in [x for x in a.genes.split(',') if x]:
        label = 'GENE:' + g
        if g in cons.index:
            r = cons.loc[g]
            rows.append(dict(marker=g, type='gene', in_graph=True, weighted_node_score_rank=int(r['rank']), nodes_of_type=len(cons),
                             in_top_fraction=bool(r.top_fraction), weighted_node_score=r.weighted_node_score,
                             exposure_alignment=r.exposure_alignment, transition_percentile=r.transition_percentile,
                             median_patient_priority_rank=nodes.median_patient_rank.get(label, np.nan),
                             cohort_priority_rank=pr.get(label, np.nan)))
        else:
            rows.append(dict(marker=g, type='gene', in_graph=False))
    for c in [x for x in a.cpgs.split(',') if x]:
        label = 'CpG:' + c
        row = dict(marker=c, type='CpG', in_graph=label in nodes.index,
                   median_patient_priority_rank=nodes.median_patient_rank.get(label, np.nan),
                   cohort_priority_rank=pr.get(label, np.nan))
        if label in cpg.index:
            r = cpg.loc[label]
            row.update(weighted_node_score_rank=int(r['rank']), nodes_of_type=len(cpg), in_top_fraction=bool(r.top_fraction),
                       weighted_node_score=r.weighted_node_score, exposure_alignment=r.exposure_alignment,
                       transition_percentile=r.transition_percentile)
        rows.append(row)
    t = pd.DataFrame(rows); t.to_csv(out / 'marker_ranks.csv', index=False); print(t.to_string(index=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', required=True); p.add_argument('--weighted', required=True); p.add_argument('--out', required=True)
    p.add_argument('--genes', default='', help='comma-separated Entrez IDs')
    p.add_argument('--cpgs', default='', help='comma-separated CpG probe IDs')
    main(p.parse_args())
