"""Exposure-learning diagnostics for the UNCHANGED ExposureGAT model.

Nothing in analyze.py, Net, the loss formula or the data loader is modified.
Net and rank_target are imported from analyze.py; the training loop below is a
line-for-line copy of analyze.fit() with extra logging and configurable weights
(defaults reproduce the current run exactly).

Part A  signal ceiling (no GNN, out-of-sample):
        repeated 5-fold ridge on (i) node-specific paired deltas, tumor and normal
        values, and (ii) node-SYMMETRIC summaries (mean/sd/quantiles across nodes).
        Net shares encoder/readout weights across all nodes and has no node
        identity, so (ii) approximates what it can reach without graph structure.
        Each is compared with permuted-label nulls. Known smoking markers are
        reported univariately when their labels exist in the input.
Part B  training diagnostics:
        per-epoch in-sample Spearman, exposure loss relative to the best constant
        prediction, per-objective gradient norms by parameter group, clipping,
        parameter displacement, and a matched permuted-label fit for every
        config/seed. Final node-rank agreement across seeds is also reported.

Stages (run in order; carry the best setting forward with --lr/--epochs):
  --stage ceiling   Part A only (minutes)
  --stage opt       learning rate x epochs at current loss weights
  --stage loss      one-at-a-time loss-weight reductions, plus exposure-only
  --stage capacity  hidden 8/16/32
All results are in-sample, exploratory diagnostics, not significance tests.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from scipy.stats import spearmanr, rankdata

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from data_loader import LungMultiOmicsDataset
from exposureGAT import Net, rank_target

GROUPS = {'encoder': ['gene', 'cpg'],
          'message_passing': ['selfmap', 'message', 'key', 'att', 'edgegate'],
          'scoring': ['nodegate', 'readout', 'bias'],
          'decoder': ['decode']}
DEFAULT_MARKERS = 'CYP1A1,1543,CYP1B1,1545,AHRR,57491,cg05575921,cg21566642,cg01940273,cg03636183'


def srho(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if np.ptp(x) == 0 or np.ptp(y) == 0: return np.nan
    return float(spearmanr(x, y).correlation)


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
def load(a):
    ds = LungMultiOmicsDataset(a.data_dir, exposure_cols=[a.exposure_col], verbose=False)
    meta = ds.meta.drop_duplicates(ds.patient_id_col).set_index(ds.patient_id_col).loc[ds.patient_ids]
    raw = meta[a.exposure_col].astype(float).to_numpy()
    keep = np.ones(len(raw), bool); reasons = []
    if a.exclude_ambiguous_zero:
        if a.status_col not in meta:
            raise ValueError(f'--exclude_ambiguous_zero needs column {a.status_col}')
        status = meta[a.status_col].astype(str).str.lower()
        never = status.str.contains('never|non-smoker|nonsmoker|non smoker|lifelong').to_numpy()
        drop = (raw == 0) & ~never
        keep &= ~drop; reasons.append(f'dropped {drop.sum()} zero-exposure patients not recorded as never-smokers')
    if a.exclude_ids:
        ids = set(Path(a.exclude_ids).read_text().split())
        drop = np.isin(ds.patient_ids, list(ids)) & keep
        keep &= ~drop; reasons.append(f'dropped {drop.sum()} listed patients')
    idx = np.where(keep)[0]
    if len(idx) < 10 or np.ptp(raw[idx]) == 0:
        raise ValueError('Need >=10 retained patients with variable exposure')
    ds.fit_scalers([ds.patient_ids[i] for i in idx])
    normals, tumors = [], []
    for i in idx:
        n, t, _, _ = ds[i]; normals.append(n.x); tumors.append(t.x)
    N = torch.stack(normals); T = torch.stack(tumors)
    info = dict(patients_total=len(raw), patients_used=int(len(idx)),
                zero_exposure_used=int((raw[idx] == 0).sum()), selection=reasons,
                genes=ds.n_genes, cpgs=ds.n_cpgs, edges=int(ds.edge_index.shape[1]))
    return ds, idx, raw[idx], N, T, info


# ----------------------------------------------------------------------------
# Part A: signal ceiling
# ----------------------------------------------------------------------------
def feature_sets(ds, N, T):
    g = ds.n_genes; n = N.numpy(); t = T.numpy(); d = t - n
    def node_specific(x):
        return np.concatenate([x[:, :g, :].reshape(len(x), -1), x[:, g:, 0]], 1)
    def symmetric(x):
        blocks = [x[:, :g, j] for j in range(3)] + ([x[:, g:, 0]] if x.shape[1] > g else [])
        out = []
        for b in blocks:
            out += [b.mean(1), b.std(1), np.abs(b).mean(1)] + list(np.percentile(b, [10, 50, 90], axis=1))
        return np.stack(out, 1)
    sets = {'node_specific_delta': node_specific(d), 'node_specific_tumor': node_specific(t),
            'node_specific_normal': node_specific(n),
            'node_symmetric_delta+tumor+normal': np.concatenate([symmetric(d), symmetric(t), symmetric(n)], 1)}
    return {k: v[:, v.std(0) > 1e-8] for k, v in sets.items()}


def cv_ridge(X, y, raw, seed, folds=5, repeats=10):
    """Mean WITHIN-test-fold Spearman over repeated K-fold.
    Pooled leave-one-out correlation is biased strongly negative (the training
    mean shifts against each held-out value), so it is not used."""
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    rng = np.random.default_rng(seed); alphas = np.logspace(-1, 6, 15); vals = []
    for _ in range(repeats):
        fold = rng.permutation(np.arange(len(y)) % folds)
        for f in range(folds):
            te = fold == f
            m = make_pipeline(StandardScaler(), RidgeCV(alphas=alphas)).fit(X[~te], y[~te])
            vals.append(srho(m.predict(X[te]), raw[te]))
    return float(np.nanmean(vals))


def ceiling(ds, N, T, raw, a, out):
    y = rank_target(raw); rng = np.random.default_rng(a.seed); rows = []
    for name, X in feature_sets(ds, N, T).items():
        obs = cv_ridge(X, y, raw, a.seed)
        null = []
        for r in range(a.ceiling_permutations):
            p = rng.permutation(len(raw)); null.append(cv_ridge(X, y[p], raw[p], a.seed + r + 1))
        null = np.array(null)
        rows.append(dict(feature_set=name, n_features=X.shape[1], cv_within_fold_spearman=obs,
                         null_mean=np.nanmean(null), null_95pct=np.nanpercentile(null, 95),
                         null_fraction_ge_observed=float(np.mean(null >= obs))))
        print(f'[ceiling] {name}: CV rho={obs:.3f}  null mean={np.nanmean(null):.3f} '
              f'95%={np.nanpercentile(null, 95):.3f}', flush=True)
    pd.DataFrame(rows).to_csv(out / 'ceiling_cv_ridge.csv', index=False)

    # Univariate descriptive associations, including known markers when present.
    g = ds.n_genes; labels = list(ds.genes) + list(ds.cpgs)
    n = N.numpy(); t = T.numpy(); recs = []
    for v, lab in enumerate(labels):
        mods = ['rna', 'meth', 'cnv'] if v < g else ['beta']
        for j, mod in enumerate(mods):
            for state, x in [('delta', t[:, v, j] - n[:, v, j]), ('tumor', t[:, v, j]), ('normal', n[:, v, j])]:
                recs.append(dict(node_label=lab, modality=mod, state=state, spearman=srho(x, raw)))
    uni = pd.DataFrame(recs).dropna()
    uni['abs_spearman'] = uni.spearman.abs()
    uni.sort_values('abs_spearman', ascending=False).head(200).to_csv(out / 'univariate_top200.csv', index=False)
    markers = {m.strip() for m in a.markers.split(',') if m.strip()}
    norm = lambda s: str(s).split('.')[0].strip()
    found = uni[uni.node_label.map(norm).isin(markers) | uni.node_label.isin(markers)]
    found.to_csv(out / 'known_marker_associations.csv', index=False)
    print(f'[ceiling] known markers found in input: {sorted(found.node_label.unique().tolist())}')


# ----------------------------------------------------------------------------
# Part B: logged training (same loop as analyze.fit)
# ----------------------------------------------------------------------------
def group_of(name):
    top = name.split('.')[0]
    return next(g for g, members in GROUPS.items() if top in members)


def train_logged(ds, N, T, raw, labels_raw, cfg, seed, a):
    torch.manual_seed(seed); np.random.seed(seed)
    m = Net(ds.n_nodes, ds.n_genes, ds.edge_index, cfg['hidden'])
    y = torch.tensor(rank_target(labels_raw), dtype=torch.float32)
    opt = torch.optim.Adam(m.parameters(), lr=cfg['lr'], weight_decay=.001)
    scale = (T - N).square().mean((0, 1)).sqrt().clamp_min(.03)
    names, params = zip(*m.named_parameters())
    init = [p.detach().clone() for p in params]
    baseline = float(F.smooth_l1_loss(torch.full_like(y, float(y.median())), y))
    log_at = set(range(1, cfg['epochs'] + 1, a.log_every)) | {cfg['epochs']} | set(a.checkpoints)
    rows = []
    for epoch in range(1, cfg['epochs'] + 1):
        m.train(); pred, c, d, recon, _, gate, edge = m(N, T)
        le = F.smooth_l1_loss(pred, y)
        err = F.smooth_l1_loss(recon / scale, (T - N) / scale, reduction='none')
        lt = (err * m.valid).sum() / (len(y) * m.valid.sum())
        terms = dict(exposure=le, transition=cfg['w_transition'] * lt,
                     node_gate=cfg['w_node_gate'] * gate.mean(), edge_gate=cfg['w_edge_gate'] * edge.mean())
        loss = sum(terms.values())
        row = None
        if epoch in log_at:
            pr = pred.detach().numpy()
            row = dict(epoch=epoch, loss=float(loss), exposure_loss=float(le), transition_loss=float(lt),
                       node_gate_mean=float(gate.mean()), edge_weight_mean=float(edge.mean()),
                       constant_baseline_loss=baseline, exposure_loss_vs_baseline=float(le) / baseline,
                       fitted_spearman_vs_supplied=srho(pr, labels_raw),
                       fitted_pearson_vs_target=float(np.corrcoef(pr, y.numpy())[0, 1]) if pr.std() > 0 else np.nan,
                       pred_sd=float(pr.std()))
            for tname, term in terms.items():
                if not term.requires_grad: continue
                grads = torch.autograd.grad(term, params, retain_graph=True, allow_unused=True)
                acc = {}
                for n_, gr in zip(names, grads):
                    if gr is not None: acc[group_of(n_)] = acc.get(group_of(n_), 0.) + float(gr.square().sum())
                for grp in GROUPS: row[f'grad_{tname}_{grp}'] = acc.get(grp, 0.) ** .5
        opt.zero_grad(); loss.backward()
        total = float(nn.utils.clip_grad_norm_(m.parameters(), 5)); opt.step()
        if row is not None:
            row['grad_total_norm_before_clip'] = total
            for grp in GROUPS:
                sel = [i for i, n_ in enumerate(names) if group_of(n_) == grp]
                num = sum(float((params[i].detach() - init[i]).square().sum()) for i in sel) ** .5
                den = sum(float(init[i].square().sum()) for i in sel) ** .5
                row[f'relative_displacement_{grp}'] = num / max(den, 1e-12)
            rows.append(row)
    m.eval()
    with torch.no_grad():
        pred, c, *_ = m(N, T)
    return pd.DataFrame(rows), c.abs().numpy(), pred.numpy()


def node_agreement(scores, k):
    """Mean across seed pairs and patients: node-rank Spearman and top-k Jaccard."""
    rh, jac = [], []
    for i in range(len(scores)):
        for j in range(i + 1, len(scores)):
            for p in range(scores[i].shape[0]):
                a_, b_ = scores[i][p], scores[j][p]
                rh.append(srho(a_, b_))
                ta = set(np.argsort(-a_)[:k]); tb = set(np.argsort(-b_)[:k])
                jac.append(len(ta & tb) / len(ta | tb))
    return (np.nanmean(rh) if rh else np.nan), (np.mean(jac) if jac else np.nan)


def configs(a):
    base = dict(hidden=a.hidden, lr=a.lr, epochs=a.epochs, w_transition=a.w_transition,
                w_node_gate=a.w_node_gate, w_edge_gate=a.w_edge_gate)
    if a.stage == 'opt':
        return {f'lr{lr:g}': {**base, 'lr': lr} for lr in a.lrs}
    if a.stage == 'loss':
        return {'current_weights': base,
                'transition_0.01': {**base, 'w_transition': .01},
                'node_gate_0.001': {**base, 'w_node_gate': .001},
                'edge_gate_0.0001': {**base, 'w_edge_gate': .0001},
                'all_three_reduced': {**base, 'w_transition': .01, 'w_node_gate': .001, 'w_edge_gate': .0001},
                'exposure_only_upper_bound': {**base, 'w_transition': 0., 'w_node_gate': 0., 'w_edge_gate': 0.}}
    if a.stage == 'capacity':
        return {f'hidden{h}': {**base, 'hidden': h} for h in a.hiddens}
    raise ValueError(a.stage)


def training(ds, N, T, raw, a, out):
    traj, finals, fits = [], [], []
    rng = np.random.default_rng(a.seed)
    perms = [rng.permutation(raw) for _ in range(a.seeds)]
    cfgs = configs(a); started = time.time(); done = 0; total = len(cfgs) * a.seeds * 2
    for cname, cfg in cfgs.items():
        scores = {'observed': [], 'permuted': []}
        for s in range(a.seeds):
            for kind, lab in [('observed', raw), ('permuted', perms[s])]:
                df, sc, pred = train_logged(ds, N, T, raw, lab, cfg, a.seed + s, a)
                df.insert(0, 'labels', kind); df.insert(0, 'seed', a.seed + s); df.insert(0, 'config', cname)
                for k, v in cfg.items(): df[k] = v
                traj.append(df); scores[kind].append(sc)
                fits.append(dict(config=cname, seed=a.seed + s, labels=kind,
                                 final_fitted_spearman=srho(pred, lab)))
                done += 1; el = time.time() - started
                print(f'[{a.stage}] {cname} seed {a.seed + s} {kind}: final rho={srho(pred, lab):.3f} '
                      f'loss/baseline={df.exposure_loss_vs_baseline.iloc[-1]:.3f}  '
                      f'({done}/{total}, ~{el / done * (total - done) / 60:.0f} min left)', flush=True)
        for kind in scores:
            r, j = node_agreement(scores[kind], min(a.top_k, ds.n_nodes))
            finals.append(dict(config=cname, labels=kind, cross_seed_node_rank_spearman=r,
                               cross_seed_topk_jaccard=j, top_k=a.top_k))
    traj = pd.concat(traj, ignore_index=True); traj.to_csv(out / 'trajectory.csv', index=False)
    pd.DataFrame(fits).to_csv(out / 'final_fits.csv', index=False)
    pd.DataFrame(finals).to_csv(out / 'node_agreement_across_seeds.csv', index=False)
    cp = traj[traj.epoch.isin(a.checkpoints + [traj.epoch.max()])]
    summ = cp.groupby(['config', 'labels', 'epoch']).agg(
        fitted_spearman_mean=('fitted_spearman_vs_supplied', 'mean'),
        fitted_spearman_sd=('fitted_spearman_vs_supplied', 'std'),
        exposure_loss_vs_baseline=('exposure_loss_vs_baseline', 'mean'),
        clip_norm=('grad_total_norm_before_clip', 'mean')).reset_index()
    wide = summ.pivot_table(index=['config', 'epoch'], columns='labels',
                            values=['fitted_spearman_mean', 'exposure_loss_vs_baseline'])
    wide.columns = ['_'.join(c) for c in wide.columns]; wide = wide.reset_index()
    if {'fitted_spearman_mean_observed', 'fitted_spearman_mean_permuted'} <= set(wide):
        wide['observed_minus_permuted_spearman'] = wide.fitted_spearman_mean_observed - wide.fitted_spearman_mean_permuted
    wide.to_csv(out / 'checkpoint_summary.csv', index=False)
    print(wide.to_string(index=False))
    plot(traj, out)


def plot(traj, out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    cfgs = list(dict.fromkeys(traj.config))
    fig, ax = plt.subplots(2, len(cfgs), figsize=(4 * len(cfgs), 7), squeeze=False)
    for j, c in enumerate(cfgs):
        for (s, kind), d in traj[traj.config == c].groupby(['seed', 'labels']):
            col = '#2886B5' if kind == 'observed' else '#A14792'
            ax[0, j].plot(d.epoch, d.fitted_spearman_vs_supplied, color=col, alpha=.7)
            ax[1, j].plot(d.epoch, d.exposure_loss_vs_baseline, color=col, alpha=.7)
        ax[0, j].set_title(c); ax[0, j].axhline(0, color='grey', lw=.5); ax[1, j].axhline(1, color='grey', lw=.5, ls='--')
        ax[1, j].set_xlabel('epoch')
    ax[0, 0].set_ylabel('in-sample Spearman\n(blue observed, purple permuted)')
    ax[1, 0].set_ylabel('exposure loss / constant baseline')
    fig.tight_layout(); fig.savefig(out / 'trajectory.png', dpi=150); plt.close(fig)

    gcols = [c for c in traj if c.startswith('grad_') and c != 'grad_total_norm_before_clip']
    g = traj[traj.labels == 'observed'].groupby(['config', 'epoch'])[gcols].mean().reset_index()
    terms = ['exposure', 'transition', 'node_gate', 'edge_gate']
    fig, ax = plt.subplots(len(GROUPS), len(cfgs), figsize=(4 * len(cfgs), 2.6 * len(GROUPS)), squeeze=False)
    for j, c in enumerate(cfgs):
        d = g[g.config == c]
        for i, grp in enumerate(GROUPS):
            for t in terms:
                col = f'grad_{t}_{grp}'
                if col in d and d[col].abs().sum() > 0: ax[i, j].semilogy(d.epoch, d[col].clip(lower=1e-12), label=t)
            ax[i, j].set_title(f'{c}: {grp}', fontsize=9)
        ax[0, j].legend(fontsize=7)
    fig.supylabel('gradient norm per objective (weighted)'); fig.tight_layout()
    fig.savefig(out / 'gradient_norms.png', dpi=150); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data_dir', required=True); p.add_argument('--out_dir', required=True)
    p.add_argument('--stage', choices=['ceiling', 'opt', 'loss', 'capacity'], required=True)
    p.add_argument('--exposure_col', default='smoking_intensity'); p.add_argument('--status_col', default='smoking_status')
    p.add_argument('--exclude_ambiguous_zero', action='store_true',
                   help='drop exposure==0 patients whose status is not never/lifelong non-smoker (incl. unreported)')
    p.add_argument('--exclude_ids', default='', help='text file of patient_ids to exclude')
    p.add_argument('--epochs', type=int, default=400); p.add_argument('--lr', type=float, default=.003)
    p.add_argument('--lrs', type=float, nargs='+', default=[.003, .001, .01])
    p.add_argument('--hidden', type=int, default=8); p.add_argument('--hiddens', type=int, nargs='+', default=[8, 16, 32])
    p.add_argument('--w_transition', type=float, default=.05); p.add_argument('--w_node_gate', type=float, default=.01)
    p.add_argument('--w_edge_gate', type=float, default=.001)
    p.add_argument('--seeds', type=int, default=3); p.add_argument('--seed', type=int, default=42)
    p.add_argument('--log_every', type=int, default=10)
    p.add_argument('--checkpoints', type=int, nargs='+', default=[60, 100, 200, 400, 800])
    p.add_argument('--top_k', type=int, default=100)
    p.add_argument('--ceiling_permutations', type=int, default=20)
    p.add_argument('--markers', default=DEFAULT_MARKERS)
    p.add_argument('--threads', type=int, default=0)
    a = p.parse_args()
    if a.threads: torch.set_num_threads(a.threads)
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()): raise ValueError('Use a new empty output directory')
    out.mkdir(parents=True, exist_ok=True)
    ds, idx, raw, N, T, info = load(a)
    info.update(config=vars(a)); (out / 'run_info.json').write_text(json.dumps(info, indent=2, default=str))
    print(json.dumps({k: v for k, v in info.items() if k != 'config'}, indent=2))
    pd.DataFrame({'patient_id': [ds.patient_ids[i] for i in idx], 'exposure': raw}).to_csv(out / 'patients_used.csv', index=False)
    if a.stage == 'ceiling': ceiling(ds, N, T, raw, a, out)
    else: training(ds, N, T, raw, a, out)


if __name__ == '__main__':
    main()
