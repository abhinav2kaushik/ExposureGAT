#!/usr/bin/env python3
"""Exposure-label controls for the weighted node score (no retraining).

Compares real-label initialization fits with shuffled-label (permutation) fits:
  1. Ranking agreement. Patient level: per patient, genes ranked by Q = W x |signed contribution|
     of each fit; cohort level: genes ranked by each fit's W. Spearman correlation and top-100
     Jaccard overlap for every real-real and real-shuffled pair of fits.
  2. Alignment counts. Number of genes whose signed contributions correlate with the fit's own
     supplied exposure at Spearman rho >= 0.2/0.3/0.4/0.5 (positive and absolute).
  3. Fixed top 5% agreement. The top 5% of genes by the real-label consensus W are fixed; for each
     permutation fit, Spearman correlation between its W and the consensus W across those genes.
Aggregate comparisons between groups of fits use two-sided run-group randomization tests
(9,999 relabelings); Holm adjustment across five primary comparisons (four agreement outcomes and
the rho >= 0.30 positive-alignment count). These compare fits, not individual genes.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


def agreement(values, pairs, idx, k):
    ranks = rankdata(values, axis=1)
    ranks -= ranks.mean(1, keepdims=True)
    norm = np.sqrt((ranks * ranks).sum(1))
    den = norm[:, None] * norm[None, :]
    corr = np.divide(ranks @ ranks.T, den, out=np.full(den.shape, np.nan), where=den > 0)
    tops = [set(idx[np.lexsort((idx, -v))[:k]]) for v in values]
    return [(i, j, float(corr[i, j]), len(tops[i] & tops[j]) / len(tops[i] | tops[j])) for i, j in pairs]


def main(a):
    model, out = Path(a.model), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    w = pd.read_csv(Path(a.weighted) / 'weighted_node_scores_by_fit.csv.gz')
    w = w[w.kind.isin(['initialization', 'permutation'])]
    runs = pd.read_csv(model / 'runs.csv')
    runs = runs[runs.kind.isin(['initialization', 'permutation'])]
    if (runs.kind == 'initialization').sum() < 2 or (runs.kind == 'permutation').sum() < 1:
        print('Label controls need >= 2 initialization and >= 1 permutation fits; skipped.')
        return
    weights, contribs, targets, names, kinds, counts, ids, idx = [], [], [], [], [], [], None, None
    for row in runs.itertuples():
        name = f'{row.kind}_{row.repeat}'
        t = w[w.run == name].sort_values('node_idx')
        if idx is not None and not np.array_equal(idx, t.node_idx.to_numpy()):
            raise ValueError('gene order mismatch')
        idx = t.node_idx.to_numpy()
        m = pd.read_csv(model / name / 'fit_membership.csv')
        if ids is None:
            ids = m.patient_id.tolist()
        if ids != m.patient_id.tolist():
            raise ValueError('patient order mismatch')
        with np.load(model / name / 'attributions.npz') as z:
            c = z['signed'][:, idx]
        rho = t.exposure_alignment.to_numpy()
        for thr in (.2, .3, .4, .5):
            for mode, mask in (('Positive alignment', rho >= thr), ('Absolute alignment', np.abs(rho) >= thr)):
                counts.append(dict(run=name, kind=row.kind, threshold=thr, mode=mode, count=int(mask.sum()), genes=len(idx)))
        weights.append(t.score.to_numpy()); contribs.append(c); targets.append(m.supplied_exposure.to_numpy(float))
        names.append(name); kinds.append(row.kind)
    ww = np.stack(weights); q = ww[:, None, :] * np.abs(np.stack(contribs))
    k = min(100, len(idx))
    real = [i for i, x in enumerate(kinds) if x == 'initialization']
    null = [i for i, x in enumerate(kinds) if x == 'permutation']
    pairs = list(itertools.combinations(real, 2)) + list(itertools.product(real, null))
    group = lambda j: 'Real vs real' if kinds[j] == 'initialization' else 'Real vs shuffled'
    pr = []
    for p, pid in enumerate(ids):
        res = agreement(q[:, p, :], pairs, idx, k)
        for g in ('Real vs real', 'Real vs shuffled'):
            s = [r for r in res if group(r[1]) == g]
            pr.append(dict(patient_id=pid, comparison=g, rank_spearman=np.nanmean([r[2] for r in s]),
                           top100_jaccard=np.mean([r[3] for r in s])))
    pd.DataFrame(pr).to_csv(out / 'patient_ranking_agreement.csv', index=False)
    pd.DataFrame([dict(run1=names[i], run2=names[j], comparison=group(j), rank_spearman=r, top100_jaccard=jc)
                  for i, j, r, jc in agreement(ww, pairs, idx, k)]).to_csv(out / 'cohort_ranking_agreement.csv', index=False)
    counts = pd.DataFrame(counts); counts.to_csv(out / 'alignment_counts.csv', index=False)

    # Run-group randomization tests.
    n = len(kinds)
    allpairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    def mats(values):
        r = np.zeros((2, n, n))
        for i, j, rr, jc in agreement(values, allpairs, idx, k):
            r[:, i, j] = r[:, j, i] = [rr, jc]
        return r
    cohort = mats(ww)
    patient = np.mean([mats(q[:, p, :]) for p in range(q.shape[1])], axis=0)
    isreal = np.array([x == 'initialization' for x in kinds])
    rng = np.random.default_rng(a.seed)
    masks = np.stack([isreal] + [rng.permutation(isreal) for _ in range(a.draws)])

    def contrast(m, g):
        rr = np.outer(g, g); np.fill_diagonal(rr, False)
        return m[rr].mean() - m[np.outer(g, ~g)].mean()
    rows = []
    for level, mm in (('patient', patient), ('cohort', cohort)):
        for metric, m in zip(('rank_spearman', 'top100_jaccard'), mm):
            stat = np.array([contrast(m, g) for g in masks])
            rows.append(dict(level=level, metric=metric, effect=stat[0],
                             p=(1 + np.sum(np.abs(stat[1:]) >= abs(stat[0]))) / (a.draws + 1)))
    for mode in counts['mode'].unique():
        for thr in sorted(counts.threshold.unique()):
            t = counts[(counts['mode'] == mode) & (counts.threshold == thr)].set_index('run')
            vals = np.array([t.loc[x, 'count'] for x in names])
            stat = np.array([vals[g].mean() - vals[~g].mean() for g in masks])
            rows.append(dict(level='counts', metric=f'{mode}: {thr}', effect=stat[0],
                             p=(1 + np.sum(np.abs(stat[1:]) >= abs(stat[0]))) / (a.draws + 1)))
    r = pd.DataFrame(rows)
    r['test'] = 'Two-sided run-group randomization'
    primary = (r.level != 'counts') | (r.metric == 'Positive alignment: 0.3')
    ix = r.index[primary]
    o = ix[np.argsort(r.loc[ix, 'p'].to_numpy())]
    adj = np.maximum.accumulate(r.loc[o, 'p'].to_numpy() * np.arange(len(o), 0, -1))
    r['p_holm'] = np.nan
    r.loc[o, 'p_holm'] = np.minimum(1, adj)
    r.to_csv(out / 'statistical_comparisons.csv', index=False)

    # Fixed top 5% genes of the real-label consensus versus each permutation fit.
    cons = pd.read_csv(Path(a.weighted) / 'weighted_node_scores.csv')
    fixed = cons.sort_values('rank').head(int(np.ceil(.05 * len(cons))))[['node_idx', 'weighted_node_score']]
    fx = []
    for i in null:
        s = pd.Series(ww[i], index=idx).loc[fixed.node_idx].to_numpy()
        rho = spearmanr(fixed.weighted_node_score, s).statistic if np.ptp(s) > 0 else np.nan
        fx.append(dict(permutation=int(names[i].split('_')[-1]), genes=len(fixed), spearman=rho))
    pd.DataFrame(fx).sort_values('permutation').to_csv(out / 'fixed_top5_correlations.csv', index=False)
    summary = dict(real_fits=len(real), permutation_fits=len(null), patients=len(ids), genes=len(idx),
                   draws=a.draws, holm_family='4 agreement outcomes + positive alignment rho >= 0.30',
                   mean_alignment_count_real=float(counts[(counts['mode'] == 'Positive alignment') & (counts.threshold == .3) & (counts.kind == 'initialization')]['count'].mean()),
                   mean_alignment_count_shuffled=float(counts[(counts['mode'] == 'Positive alignment') & (counts.threshold == .3) & (counts.kind == 'permutation')]['count'].mean()),
                   fixed_top5_mean_spearman=float(np.nanmean([x['spearman'] for x in fx])))
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(r.to_string(index=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', required=True); p.add_argument('--weighted', required=True); p.add_argument('--out', required=True)
    p.add_argument('--draws', type=int, default=9999); p.add_argument('--seed', type=int, default=20261103)
    main(p.parse_args())
