#!/usr/bin/env python3
"""Summary tables for the simulation benchmark (independent simulated datasets): AUROCs of Set A
versus Set C (and secondary comparisons), top-k recovery, single-run AUROCs under initialization,
patient-subset and exposure-permutation fits, and design checks (latent magnitude invariance,
equal node degree). Figures are drawn by src/figures/simulation_figures.R.
Ground truth is used only for evaluation. One dataset seed is the unit of replication.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from simulation_scaling import SimulationDataset

ORDER=['A','C','B']
PAIRS=[('A','C'),('A','B'),('C','B')]
def rho(x,y):
    return float(spearmanr(x,y).statistic) if np.ptp(x)>1e-10 and np.ptp(y)>1e-10 else 0.
def auc(labels,score,a,b):
    mask=np.isin(labels,[a,b]); y=labels[mask]==a;s=np.asarray(score)[mask]
    assert np.isfinite(s).all() and y.any() and (~y).any()
    return float((rankdata(s)[y].sum()-y.sum()*(y.sum()+1)/2)/(y.sum()*(~y).sum()))
def ranks(x):return np.apply_along_axis(lambda v:rankdata(-v),-1,x)
def main(a):
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(a.manifest)
    if manifest.seed.duplicated().any():raise ValueError('Use one independent dataset per seed')
    metrics=[];tops=[];design=[];runrows=[];nodeframes=[];patientframes=[];degreeframes=[]
    example=None
    for row in manifest.itertuples():
        data=Path(row.data_dir);run=Path(row.analysis_dir)
        gt=pd.read_csv(data/'sim_ground_truth_nodes.csv').sort_values('node_idx')
        conf=json.loads((run/'config.json').read_text())
        if conf['scaling']!='shared':raise ValueError('This benchmark requires shared modality scaling')
        cohort=pd.read_csv(run/'cohort_summary.csv').merge(gt,on=['node_idx','node_label'],validate='one_to_one').sort_values('node_idx')
        if len(cohort)!=len(gt):raise ValueError('Incomplete truth/score join')
        labels=cohort.set_membership.to_numpy();v=len(cohort)
        priority=-cohort.median_patient_rank.to_numpy()
        cohort['normalized_priority']=(v-cohort.median_patient_rank)/(v-1)
        cohort['seed']=row.seed;nodeframes.append(cohort)
        for aa,bb in PAIRS:metrics.append(dict(seed=row.seed,comparison=f'{aa} vs {bb}',AUROC=auc(labels,priority,aa,bb)))
        k=min(a.top_k,v)
        chosen=cohort.sort_values(['median_patient_rank','node_idx']).head(k)
        for c in ORDER:
            count=int((chosen.set_membership==c).sum());size=int((labels==c).sum())
            tops.append(dict(seed=row.seed,set=c,k=k,selected=count,set_size=size,recall=count/size,expected=k*size/v,fold_enrichment=count/(k*size/v)))
        cohort.sort_values(['median_patient_rank','node_idx']).head(min(100,v)).to_csv(out/f'top100_seed_{row.seed}.csv',index=False)
        ds=SimulationDataset(str(data),exposure_cols=[conf['exposure_col']],verbose=False)
        ds.fit_scalers(ds.patient_ids)
        exposed=np.array([float(ds.meta[ds.meta.patient_id==p][conf['exposure_col']].iloc[0]) for p in ds.patient_ids])
        membership=pd.read_csv(run/'initialization_1/fit_membership.csv')
        assert list(membership.patient_id)==list(ds.patient_ids)
        delta=np.stack([(ds[i][1].x-ds[i][0].x).numpy() for i in range(len(ds))])
        latent=np.load(data/'sim_latent_node_changes.npz')['delta']
        l1range=np.ptp(np.abs(latent).sum(1));l2range=np.ptp((latent**2).sum(1))
        assert max(l1range,l2range)<1e-8,'Latent magnitude invariance failed'
        edges=pd.read_csv(data/'edges_ppi.csv');und=np.unique(np.sort(edges[['from_idx','to_idx']].to_numpy(),axis=1),axis=0)
        und=und[und[:,0]!=und[:,1]]
        degrees=np.bincount(und.ravel(),minlength=v)
        assert np.ptp(degrees)==0,'Equal-degree design check failed'
        degreeframes.append(pd.DataFrame(dict(seed=row.seed,node_idx=np.arange(v),set=labels,degree=degrees)))
        touches=np.any(labels[und]=='A',axis=1).mean()
        norms=pd.DataFrame(dict(seed=row.seed,patient_id=ds.patient_ids,exposure=exposed,transformed_L1=np.abs(delta).sum((1,2)),transformed_L2_squared=(delta**2).sum((1,2))))
        for m in ['rna','meth','cnv']:
            x=pd.read_csv(data/f'X_{m}_gene_features.csv',index_col=0)
            d=x.loc[[p+'-TP' for p in ds.patient_ids]].to_numpy()-x.loc[[p+'-NT' for p in ds.patient_ids]].to_numpy()
            norms[m+'_L1']=np.abs(d).sum(1);norms[m+'_L2_squared']=(d*d).sum(1)
        for name in norms.columns[3:]:design.append(dict(seed=row.seed,measure=name,exposure_rho=rho(norms[name],exposed),latent_L1_range=l1range,latent_L2_squared_range=l2range,degree_min=degrees.min(),degree_max=degrees.max(),fraction_all_edges_touching_A=touches))
        patientframes.append(norms)
        for rr in pd.read_csv(run/'runs.csv').itertuples():
            z=np.load(run/f'{rr.kind}_{rr.repeat}'/'attributions.npz')
            r=ranks(z['score']);s=-np.median(r,axis=0)
            runrows.append(dict(seed=row.seed,kind=rr.kind,repeat=rr.repeat,A_vs_C_AUROC=auc(labels,s,'A','C')))
        if example is None:example=row
    met=pd.DataFrame(metrics);top=pd.DataFrame(tops);audit=pd.DataFrame(design);runs=pd.DataFrame(runrows)
    nodes=pd.concat(nodeframes);patients=pd.concat(patientframes);deg=pd.concat(degreeframes)
    for name,frame in [('AUROC_by_dataset',met),('topk_recovery',top),('design_audit',audit),('repetition_AUROC',runs),('node_scores',nodes),('patient_magnitude_audit',patients),('actual_node_degrees',deg)]:frame.to_csv(out/(name+'.csv'),index=False)
    summary=met.groupby('comparison').AUROC.agg(['count','mean','std','min','max'])
    summary.to_csv(out/'AUROC_summary.csv')
    row=example
    # Single-run A-versus-C AUROC: average runs within each dataset, then mean and SD across datasets.
    rr=runs.groupby(['seed','kind'],as_index=False).A_vs_C_AUROC.mean()
    rr.groupby('kind').A_vs_C_AUROC.agg(['count','mean','std','min','max']).to_csv(out/'repetition_AUROC_summary.csv')
    top.groupby('set')[['selected','expected','recall','fold_enrichment']].mean().to_csv(out/'topk_recovery_summary.csv')
    (out/'analysis_notes.json').write_text(json.dumps(dict(independent_datasets=len(manifest),representative_dataset=int(row.seed),score='negative median_patient_rank',top_k=a.top_k,uncertainty='Across independent datasets: SD, not node-wise standard error',edge_recovery_evaluated=False,cpg_relationships_evaluated=False,warning='In-sample patient score correlation is an optimization diagnostic, not evidence of recovery. No automatic pass/fail or p-values.'),indent=2))
    print(summary.to_string())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);p.add_argument('--top-k',type=int,default=50)
    main(p.parse_args())
