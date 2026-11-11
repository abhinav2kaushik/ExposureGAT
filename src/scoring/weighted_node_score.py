#!/usr/bin/env python3
"""Weighted node score: the main exposure-aligned node ranking of ExposureGAT.

Computed for every node type in the graph (gene nodes, CpG nodes, ...); each type is scored and
ranked separately. For every fitted model (initialization, patient-subset and permutation fits) and
every node:
  A = exposure alignment: Spearman correlation across patients between the gene's signed
      contribution c(p,i) and the exposure supplied to that fit (recorded exposure for real-label
      fits, the shuffled values for permutation fits).
  M = paired-transition magnitude: median across patients of ||h_tumor - h_normal||_2.
  T = percentile of M among nodes of the same type, (rank(M) - 0.5) / N.
  W = max(A, 0) * ((1 - gamma) + gamma * T), gamma = 0.10 by default.
The cohort (consensus) score takes the median A and median M across the real-label initialization
fits and then recomputes T and W. Nodes are ranked by W within type (ties by node index).

Patient-level view: Q(p,i) = W(i) * |median signed contribution of node i in patient p across
real-label initialization fits|. For each patient, the top fraction of nodes by Q is exported.

Outputs (in --out), gene nodes: weighted_node_scores.csv (main output), weighted_node_scores_by_fit.csv.gz,
patient_weighted_contributions.csv.gz, top_genes_signed_contributions.csv, top_gene_alignment_bins.csv
and two figures (top-gene heatmap; exposure alignment of the top 5/10/20% of genes). Other node types
get the same files with the type as suffix (e.g. weighted_node_scores_CpG.csv,
top_CpGs_signed_contributions.csv). settings.json lists all node types.
No retraining; reads the saved fits written by exposuregat_model.py.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def alignment(contrib, exposure):
    """Spearman correlation of each column of contrib (patients x genes) with exposure."""
    r = rankdata(contrib, axis=0)
    r -= r.mean(0)
    t = rankdata(exposure)
    t -= t.mean()
    den = np.sqrt((r * r).sum(0) * (t * t).sum())
    return np.divide((r * t[:, None]).sum(0), den, out=np.zeros(contrib.shape[1]), where=den > 0)


def weighted(a, m, gamma):
    t = (rankdata(m, method='average') - .5) / len(m)
    return np.maximum(a, 0) * (1 - gamma + gamma * t), t


def order(score, idx):
    return np.lexsort((idx, -score))


def load_fits(model, node_idx):
    runs = pd.read_csv(model / 'runs.csv')
    fits, ids = [], None
    for row in runs.itertuples():
        name = f'{row.kind}_{row.repeat}'
        folder = model / name
        mem = pd.read_csv(folder / 'fit_membership.csv')
        if ids is None:
            ids = mem.patient_id.tolist()
        if mem.patient_id.tolist() != ids:
            raise ValueError(f'patient order differs in {name}')
        y = mem.supplied_exposure.to_numpy(float)
        if not np.isfinite(y).all() or np.ptp(y) == 0:
            raise ValueError(f'invalid exposure in {name}')
        with np.load(folder / 'attributions.npz') as z:
            c = z['signed'][:, node_idx]
        m = np.load(folder / 'node_transition.npy', mmap_mode='r')[:, node_idx]
        if not (np.isfinite(c).all() and np.isfinite(m).all()):
            raise ValueError(f'non-finite values in {name}')
        fits.append(dict(kind=row.kind, run=name, exposure=y, signed=np.asarray(c),
                         transition=np.median(m, axis=0)))
    return fits, ids


def plot_heatmap(matrix, genes, exposure, scores, path, noun='genes'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sd = matrix.std(1, keepdims=True)
    z = (matrix - matrix.mean(1, keepdims=True)) / np.where(sd > 0, sd, 1)
    n = len(genes)
    fig = plt.figure(figsize=(13, max(6, .24 * n + 2)))
    gs = fig.add_gridspec(2, 3, height_ratios=[.5, .24 * n], width_ratios=[10, .25, 1.6], hspace=.08, wspace=.15)
    ax0 = fig.add_subplot(gs[0, 0])
    im0 = ax0.imshow(exposure[None, :], aspect='auto', cmap='Greys')
    ax0.set_xticks([]); ax0.set_yticks([0], ['Recorded exposure'])
    fig.colorbar(im0, cax=fig.add_subplot(gs[0, 1]))
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(z, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)
    ax.set_yticks(np.arange(n), genes, fontsize=8); ax.set_xticks([])
    ax.set_xlabel('Patients ordered by recorded exposure (low to high)')
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 1]), extend='both')
    cb.set_label('Signed contribution (row z-score)')
    bx = fig.add_subplot(gs[1, 2])
    bx.barh(np.arange(n), scores, color='#33869D'); bx.set_ylim(n - .5, -.5); bx.set_yticks([])
    bx.set_xlabel('Weighted\nnode score')
    fig.suptitle(f'Patient-specific contributions of the top {n} {noun} by weighted node score', fontweight='bold')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight'); fig.savefig(path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_alignment(table, path, noun='genes'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {'Top 5%': '#CF5D3A', 'Top 10%': '#287D8E', 'Top 20%': '#8263A7'}
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, d in table.groupby('selection', sort=False):
        ax.plot(d.bin_mid, d.percent, marker='o', label=label, color=colors.get(label))
    ax.set_xlim(0, 1); ax.set_ylim(0, 100)
    ax.set_xlabel('Absolute Spearman correlation (0.1-wide bins)'); ax.set_ylabel(f'Selected {noun} per bin (%)')
    ax.set_title(f'Exposure alignment of top-ranked {noun}', fontweight='bold'); ax.legend(frameon=False); ax.grid(alpha=.2)
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight'); fig.savefig(path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def score_type(a, out, nodes, fits, ids, idcol, suffix, noun):
    """Weighted node score for one node type; fits hold that type's columns only."""
    idx = nodes.node_idx.to_numpy()
    k = int(np.ceil(a.top_fraction * len(nodes)))
    real = [f for f in fits if f['kind'] == 'initialization']
    rows = []
    for f in fits:
        f['alignment'] = alignment(f['signed'], f['exposure'])
        w, t = weighted(f['alignment'], f['transition'], a.gamma)
        sel = np.zeros(len(idx), bool); sel[order(w, idx)[:k]] = True
        f['score'] = w
        rows.append(nodes.assign(kind=f['kind'], run=f['run'], score=w, exposure_alignment=f['alignment'],
                                 transition_magnitude=f['transition'], transition_percentile=t, selected=sel))
    pd.concat(rows, ignore_index=True).to_csv(out / f'weighted_node_scores_by_fit{suffix}.csv.gz', index=False)

    a_med = np.median([f['alignment'] for f in real], axis=0)
    m_med = np.median([f['transition'] for f in real], axis=0)
    w, t = weighted(a_med, m_med, a.gamma)
    o = order(w, idx)
    cons = nodes.iloc[o].copy()
    cons.insert(0, 'rank', np.arange(1, len(o) + 1))
    cons['weighted_node_score'] = w[o]; cons['exposure_alignment'] = a_med[o]
    cons['transition_magnitude'] = m_med[o]; cons['transition_percentile'] = t[o]
    cons['top_fraction'] = cons['rank'] <= k
    cons.to_csv(out / f'weighted_node_scores{suffix}.csv', index=False)

    # Patient-level weighted contributions.
    c_med = np.median([f['signed'] for f in real], axis=0)
    q = np.abs(c_med) * w[None, :]
    prow = []
    for i, pid in enumerate(ids):
        for r, j in enumerate(order(q[i], idx)[:k], 1):
            prow.append(dict(patient_id=pid, rank=r, node_idx=int(idx[j]), node_label=nodes.node_label.iloc[j],
                             weighted_node_score=w[j], signed_contribution=c_med[i, j], weighted_contribution=q[i, j]))
    pd.DataFrame(prow).to_csv(out / f'patient_weighted_contributions{suffix}.csv.gz', index=False)

    # Top-node heatmap data (signed median contributions; patients ordered by recorded exposure).
    exposure = real[0]['exposure']
    po = np.argsort(exposure, kind='stable')
    top = o[:min(a.heatmap_genes, len(o))]
    labels = nodes[idcol].to_numpy()[top]
    pd.DataFrame(c_med[po][:, top].T, index=nodes.node_label.to_numpy()[top],
                 columns=np.array(ids)[po]).rename_axis('node_label').to_csv(out / f'top_{noun}_signed_contributions.csv')

    # Exposure alignment of the top 5/10/20% nodes in each real-label fit, 0.1-wide bins.
    bins = np.linspace(0, 1, 11)
    hist = []
    for frac in (.05, .10, .20):
        n = int(np.ceil(frac * len(idx)))
        per = [np.histogram(np.abs(f['alignment'][order(f['score'], idx)[:n]]), bins=bins)[0] / n * 100 for f in real]
        share = np.mean([(np.abs(f['alignment'][order(f['score'], idx)[:n]]) >= .3).mean() for f in real])
        for b, v in enumerate(np.mean(per, axis=0)):
            hist.append(dict(selection=f'Top {int(frac*100)}%', **{noun: n}, bin_low=bins[b], bin_mid=bins[b] + .05,
                             percent=v, mean_share_abs_rho_ge_0_30=share))
    hist = pd.DataFrame(hist)
    stem = 'gene' if noun == 'genes' else noun.rstrip('s')
    hist.to_csv(out / f'top_{stem}_alignment_bins.csv', index=False)
    if not a.no_plots:
        title = 'Top_genes' if noun == 'genes' else f'Top_{noun}'
        plot_heatmap(c_med[po][:, top].T, labels, exposure[po], w[top], out / f'{title}_signed_contributions_heatmap', noun)
        plot_alignment(hist, out / f'Top_{stem}_exposure_alignment', noun)
    print(f'Weighted node score ({noun}): {len(idx)} nodes, top {k} exported.')
    return dict(nodes=len(idx), top_k=k)


def main(a):
    model, out = Path(a.model), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    nodes = pd.read_csv(model / 'cohort_summary.csv').sort_values('node_idx')[['node_idx', 'node_label']]
    nodes['prefix'] = nodes.node_label.str.split(':', n=1).str[0]
    all_idx = nodes.node_idx.to_numpy()
    fits, ids = load_fits(model, all_idx)
    real = [f for f in fits if f['kind'] == 'initialization']
    if not real:
        raise ValueError('no initialization fits found')
    pos = {n: i for i, n in enumerate(all_idx)}
    prefixes = sorted(nodes.prefix.unique(), key=lambda x: (x != 'GENE', x))
    types = {}
    for prefix in prefixes:
        sub = nodes[nodes.prefix == prefix].reset_index(drop=True)
        cols = np.array([pos[n] for n in sub.node_idx])
        gene = prefix == 'GENE'
        idcol = 'entrez_id' if gene else 'feature_id'
        sub = sub.drop(columns='prefix').assign(**{idcol: sub.node_label.str.split(':', n=1).str[1]})
        noun = 'genes' if gene else f'{prefix}s'
        type_fits = [dict(kind=f['kind'], run=f['run'], exposure=f['exposure'], signed=f['signed'][:, cols],
                          transition=f['transition'][cols]) for f in fits]
        types[prefix] = score_type(a, out, sub, type_fits, ids, idcol, '' if gene else f'_{prefix}', noun)

    (out / 'settings.json').write_text(json.dumps(dict(
        gamma=a.gamma, top_fraction=a.top_fraction, node_types=types, real_label_fits=len(real),
        fits=len(fits), formula='W = max(A,0) * ((1-gamma) + gamma*T); scored and ranked within each node type',
        consensus='median A and median M across real-label initialization fits, then T and W recomputed',
        patient_view='Q = W * |median signed contribution across real-label initialization fits|'), indent=2))
    print(f'{len(real)} real-label and {len(fits) - len(real)} other fits.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', required=True, help='model output folder (results/<name>/model)')
    p.add_argument('--out', required=True)
    p.add_argument('--gamma', type=float, default=.10)
    p.add_argument('--top_fraction', type=float, default=.10)
    p.add_argument('--heatmap_genes', type=int, default=30, help='nodes shown in each top-node heatmap')
    p.add_argument('--no-plots', dest='no_plots', action='store_true')
    a = p.parse_args()
    if not (0 <= a.gamma <= 1 and 0 < a.top_fraction <= 1):
        p.error('gamma must be in [0,1] and top_fraction in (0,1]')
    main(a)
