#!/usr/bin/env python3
"""Representation of predefined MSigDB gene signatures among genes ranked by the weighted node score.

For each fit (real-label initialization, permutation and patient-subset fits), the top fraction of
gene nodes by that fit's weighted node score (default top 10%) is intersected with each retained
signature of the catalog (resources/signature_catalog.csv; column included = True). Reported per
signature and fit type: mean and SD of selected genes, % of the signature's genes present in the
model, and fold relative to the selection fraction. Signature genes are also pooled into unique
lung cancer-related and smoking-related gene sets (signatures labelled "Smoking and lung cancer"
count in both) and their representation in the top 5/10/20% of real-label fits is summarized.
Descriptive comparisons; no enrichment P values.

Inputs: weighted node score folder (weighted_node_scores_by_fit.csv.gz), a folder of MSigDB GMT
files with Entrez IDs (msigdb.v2026.1.Hs.entrez.gmt, plain or zipped; the repository ships it as
MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip) and the
signature catalog.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

COLLECTION = [('GOBP_', 'C5 (GO:BP)'), ('HP_', 'C5 (HPO)'), ('KEGG_MEDICUS', 'C2 (CP:KEGG_MEDICUS)'),
              ('KEGG_', 'C2 (CP:KEGG_LEGACY)'), ('WP_', 'C2 (CP:WIKIPATHWAYS)'), ('REACTOME_', 'C2 (CP:REACTOME)')]


def collection(name):
    return next((c for p, c in COLLECTION if name.startswith(p)), 'C2 (CGP)')


def gmt_lines(gmt_dir):
    """Lines of every GMT file in gmt_dir: plain .gmt, or zipped (.gmt.zip / .zip) as stored in the repository."""
    import zipfile
    files = sorted(Path(gmt_dir).glob('*.gmt')) + sorted(Path(gmt_dir).glob('*.zip'))
    if not files:
        raise FileNotFoundError(f'no .gmt or .gmt.zip files in {gmt_dir}')
    for f in files:
        if f.suffix == '.zip':
            with zipfile.ZipFile(f) as z:
                for name in z.namelist():
                    if name.endswith('.gmt') and not name.startswith('__MACOSX'):
                        with z.open(name) as handle:
                            for raw in handle:
                                yield raw.decode('utf-8')
        else:
            with f.open() as handle:
                yield from handle


def read_sets(gmt_dir, wanted):
    sets = {}
    for line in gmt_lines(gmt_dir):
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) < 3 or parts[0] not in wanted:
            continue
        genes = set(parts[2:])
        if parts[0] in sets and sets[parts[0]] != genes:
            raise ValueError(f'conflicting GMT definitions for {parts[0]}')
        sets[parts[0]] = genes
    missing = wanted - set(sets)
    if missing:
        raise ValueError(f'{len(missing)} catalog signatures not found in the GMT files, e.g. {sorted(missing)[:3]}')
    return sets


def main(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(a.catalog)
    cat = cat[cat.included.astype(str).str.lower() == 'true'].copy()
    sets = read_sets(a.gmt, set(cat.pathway))
    fits = pd.read_csv(Path(a.weighted) / 'weighted_node_scores_by_fit.csv.gz', dtype={'entrez_id': str})
    genes = fits[['node_idx', 'entrez_id']].drop_duplicates().sort_values('node_idx')
    n_genes = len(genes)
    k = int(np.ceil(a.top_fraction * n_genes)); background = k / n_genes
    model_genes = set(genes.entrez_id)

    rows = []
    for run, d in fits.groupby('run', sort=False):
        d = d.sort_values(['score', 'node_idx'], ascending=[False, True])
        chosen = set(d.entrez_id.iloc[:k])
        for r in cat.itertuples():
            present = sets[r.pathway] & model_genes
            n = len(present & chosen)
            rows.append(dict(kind=d.kind.iloc[0], run=run, pathway=r.pathway, category=r.category,
                             scored_genes=len(present), selected=n,
                             fraction=n / len(present) if present else np.nan, background=background))
    counts = pd.DataFrame(rows); counts.to_csv(out / 'signature_counts_by_fit.csv', index=False)
    summary = counts.groupby(['kind', 'pathway', 'category', 'scored_genes'], dropna=False).agg(
        mean_selected=('selected', 'mean'), sd_selected=('selected', 'std'), mean_fraction=('fraction', 'mean'),
        sd_fraction=('fraction', 'std'), runs=('run', 'nunique')).reset_index()
    summary['mean_fold'] = summary.mean_fraction / background
    summary.to_csv(out / 'signature_summary.csv', index=False)
    real = summary[summary.kind == 'initialization'].drop(columns='kind')
    null = summary[summary.kind == 'permutation'].drop(columns='kind')
    if len(null):
        diff = real.merge(null, on=['pathway', 'category', 'scored_genes'], suffixes=('_real', '_shuffled'))
        diff['fraction_difference'] = diff.mean_fraction_real - diff.mean_fraction_shuffled
        diff.to_csv(out / 'real_minus_shuffled.csv', index=False)

    # Supplementary table: one row per signature.
    wide = summary.pivot_table(index='pathway', columns='kind', values=['mean_selected', 'sd_selected', 'mean_fraction'])
    wide.columns = [f'{v}_{k}' for v, k in wide.columns]
    t = cat[['pathway', 'category', 'gmt_genes', 'scored_genes', 'source_url']].merge(wide.reset_index(), on='pathway', how='left')
    table = pd.DataFrame({'Gene set': t.pathway, 'MSigDB collection': t.pathway.map(collection), 'Category': t.category,
                          'Genes in gene set': t.gmt_genes, 'Genes in model': t.scored_genes})
    for kind, label in (('initialization', 'Real labels'), ('permutation', 'Shuffled labels'), ('subset', 'Patient subsets')):
        if f'mean_selected_{kind}' in t:
            table[f'{label}: mean selected'] = t[f'mean_selected_{kind}'].round(2)
            table[f'{label}: SD'] = t[f'sd_selected_{kind}'].round(2)
            table[f'{label}: % selected'] = (100 * t[f'mean_fraction_{kind}']).round(1)
    if 'mean_fraction_permutation' in t:
        table['Difference in % selected (real − shuffled)'] = (100 * (t.mean_fraction_initialization - t.mean_fraction_permutation)).round(1)
    table['Source'] = t.source_url
    table.sort_values(['Category', 'Gene set']).to_csv(out / 'gene_signature_table.csv', index=False)

    # Unique lung cancer-related and smoking-related gene sets.
    links = []
    for r in cat.itertuples():
        cats = ['Lung cancer-related'] if r.category == 'Lung-cancer-related' else ['Smoking-related']
        if r.category == 'Smoking and lung cancer':
            cats.append('Lung cancer-related')
        for c in cats:
            links += [dict(category=c, pathway=r.pathway, entrez_id=g) for g in sets[r.pathway]]
    links = pd.DataFrame(links).drop_duplicates()
    cons_file = Path(a.weighted) / 'weighted_node_scores.csv'
    if cons_file.exists():
        cons = pd.read_csv(cons_file, dtype={'entrez_id': str})[['entrez_id', 'rank', 'weighted_node_score', 'exposure_alignment', 'top_fraction']]
        m = links.merge(cons.rename(columns={'rank': 'weighted_node_score_rank', 'top_fraction': 'in_top_fraction'}),
                        on='entrez_id', how='left')
        m['weighted_node_score_rank'] = m.weighted_node_score_rank.astype('Int64')
        m.sort_values(['category', 'pathway', 'weighted_node_score_rank', 'entrez_id'], na_position='last') \
            .to_csv(out / 'gene_set_membership.csv', index=False)
    else:
        links.to_csv(out / 'gene_set_membership.csv', index=False)
    union = []
    realfits = fits[fits.kind == 'initialization']
    for c, g in links.groupby('category'):
        members = set(g.entrez_id); in_model = members & model_genes
        for frac in (.05, .10, .20):
            n = int(np.ceil(frac * n_genes))
            sel = [len(set(d.sort_values(['score', 'node_idx'], ascending=[False, True]).entrez_id.iloc[:n]) & in_model)
                   for _, d in realfits.groupby('run')]
            union.append(dict(category=c, unique_genes=len(members), genes_in_model=len(in_model), selection=f'Top {int(frac*100)}%',
                              selected_genes_per_fit=n, mean_selected=np.mean(sel), sd_selected=np.std(sel, ddof=1),
                              expected_random=len(in_model) * n / n_genes, fold=np.mean(sel) / (len(in_model) * n / n_genes)))
    pd.DataFrame(union).to_csv(out / 'gene_set_union_summary.csv', index=False)
    print(pd.DataFrame(union).round(2).to_string(index=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--weighted', required=True); p.add_argument('--gmt', required=True)
    p.add_argument('--catalog', default=str(Path(__file__).resolve().parents[2] / 'resources' / 'signature_catalog.csv'))
    p.add_argument('--out', required=True); p.add_argument('--top_fraction', type=float, default=.10)
    main(p.parse_args())
