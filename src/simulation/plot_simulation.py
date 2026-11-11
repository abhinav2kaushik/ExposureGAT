#!/usr/bin/env python3
"""Figure 2 and Supplementary Figures 1-5 data for the simulation benchmark (ten independent datasets).
Audit and plot personalized GNN simulations. No ground truth is used for training.
One dataset seed is the unit of replication; optimizer runs are NOT independent datasets.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from simulation_scaling import SimulationDataset

COL={'A':'#ed7658','B':'#939d9d','C':'#8361bb'}
ORDER=['A','C','B']
PREVIEW=False
PAIRS=[('A','C'),('A','B'),('C','B')]
def rho(x,y):
    return float(spearmanr(x,y).statistic) if np.ptp(x)>1e-10 and np.ptp(y)>1e-10 else 0.
def auc(labels,score,a,b):
    mask=np.isin(labels,[a,b]); y=labels[mask]==a;s=np.asarray(score)[mask]
    assert np.isfinite(s).all() and y.any() and (~y).any()
    return float((rankdata(s)[y].sum()-y.sum()*(y.sum()+1)/2)/(y.sum()*(~y).sum()))
def ranks(x):return np.apply_along_axis(lambda v:rankdata(-v),-1,x)
def save(fig,out,name):
    if PREVIEW:fig.text(.5,.005,'SHORT-RUN SOFTWARE CHECK — NOT MANUSCRIPT RESULTS',ha='center',color='#a03030',fontsize=9)
    fig.savefig(out/(name+'.png'),dpi=200,bbox_inches='tight')
    fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def tidy(ax,title,x=None,y=None):
    ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
    if x:ax.set_xlabel(x)
    if y:ax.set_ylabel(y)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.18);ax.set_axisbelow(True)
def dots(ax,frame,category,value,categories):
    for j,c in enumerate(categories):
        v=frame.loc[frame[category]==c,value].to_numpy()
        ax.scatter(j+np.linspace(-.10,.10,len(v)),v,s=26,alpha=.75)
        if len(v):ax.errorbar(j,np.mean(v),yerr=np.std(v,ddof=1) if len(v)>1 else 0,fmt='k_',capsize=4)
    ax.set_xticks(range(len(categories)),categories)

def main(a):
    global PREVIEW
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(a.manifest)
    if manifest.seed.duplicated().any():raise ValueError('Use one independent dataset per seed')
    metrics=[];tops=[];design=[];runrows=[];nodeframes=[];patientframes=[];degreeframes=[]
    example=None
    for row in manifest.itertuples():
        data=Path(row.data_dir);run=Path(row.analysis_dir)
        gt=pd.read_csv(data/'sim_ground_truth_nodes.csv').sort_values('node_idx')
        conf=json.loads((run/'config.json').read_text())
        PREVIEW=PREVIEW or conf['epochs']<5
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
        if example is None:
            scores=np.load(run/'initialization_1/attributions.npz')
            pn=pd.read_csv(run/'patient_node_scores.csv')
            example=(row,data,run,cohort,norms,scores,pn,conf)
    met=pd.DataFrame(metrics);top=pd.DataFrame(tops);audit=pd.DataFrame(design);runs=pd.DataFrame(runrows)
    nodes=pd.concat(nodeframes);patients=pd.concat(patientframes);deg=pd.concat(degreeframes)
    for name,frame in [('AUROC_by_dataset',met),('topk_recovery',top),('design_audit',audit),('repetition_AUROC',runs),('node_scores',nodes),('patient_magnitude_audit',patients),('actual_node_degrees',deg)]:frame.to_csv(out/(name+'.csv'),index=False)
    summary=met.groupby('comparison').AUROC.agg(['count','mean','std','min','max'])
    summary.to_csv(out/'AUROC_summary.csv')
    row,data,run,cohort,norms,z,pn,conf=example
    cfg=json.loads((data/'simulation_config.json').read_text());n=len(cohort)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,3,height_ratios=[.8,1],hspace=.5,wspace=.35)
    ax=fig.add_subplot(gs[0,:]);ax.axis('off');ax.set_title('A. Paired multi-omics simulation and exposure-aware analysis',loc='left',fontweight='bold')
    texts=[f"{cfg['n_patients']} paired subjects / dataset\n{cfg['n_genes']} gene nodes\nRNA · methylation · CNV\nMatched, fixed graph topology\nDegree = {cfg['degree']} for every node",
           f"Set A (n={cfg['n_a']}): exposure-linked\nredistribution of paired changes\nSet C (n={cfg['n_c']}): fixed transition changes\nSet B: background variation\nLatent total L1 and L2² are constant",
           'Attention and gated message passing\nExposure-rank regression\nAuxiliary transition reconstruction\nPatient-specific contextual contributions\nInitialization · subsets · permuted exposure\nGround truth used only for evaluation']
    for i,t in enumerate(texts):
        ax.text(.015+i*.335,.50,t,transform=ax.transAxes,va='center',fontsize=9,bbox=dict(boxstyle='round,pad=.8',fc='#f4f7fa',ec='#687984'))
    ax=fig.add_subplot(gs[1,0]);e=rankdata(norms.exposure)/len(norms)
    ax.scatter(e,z['patient_score'],s=16,c='#323c48',alpha=.7)
    tidy(ax,'B. Fitted patient exposure score','Exposure percentile (ties averaged)','Patient score (standardized rank units)')
    ax.text(.03,.95,f"Descriptive, in-sample\nSpearman ρ = {rho(e,z['patient_score']):.2f}\nDataset seed {row.seed}; initialization 1",transform=ax.transAxes,va='top')
    ax=fig.add_subplot(gs[1,1]);dots(ax,met,'comparison','AUROC',['A vs C','A vs B','C vs B'])
    ax.axhline(.5,ls='--',c='gray');ax.set_ylim(0,1.03);tidy(ax,'C. Cohort node recovery',y='AUROC (higher priority = positive set)')
    ax.text(.02,.02,'Dots: independent datasets; bars: mean ± SD',transform=ax.transAxes,fontsize=8)
    ax=fig.add_subplot(gs[1,2])
    for j,c in enumerate(ORDER):
        f=top[top['set']==c];ax.bar(j,f.selected.mean(),color=COL[c],alpha=.7)
        ax.scatter(j+np.linspace(-.1,.1,len(f)),f.selected,c='black',s=16)
        ax.scatter(j,f.expected.mean(),marker='_',s=400,c='black',label='Random-selection expectation' if j==0 else None)
    ax.set_xticks(range(3),['Set '+x for x in ORDER]);tidy(ax,f'D. Top-{min(a.top_k,n)} cohort-ranked nodes',y='Number selected');ax.legend(fontsize=7)
    fig.suptitle(f"ExposureGAT simulation evaluation — {len(manifest)} independent datasets",fontsize=14)
    save(fig,out,'Figure_2')
    # Training histories from all repetitions of the first predetermined dataset.
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    for kind,color in [('initialization','#2886b5'),('subset','#c17a26'),('permutation','#a14792')]:
        hs=[pd.read_csv(p) for p in sorted(run.glob(kind+'_*/history.csv'))]
        for ax,key in zip(axes,['loss','exposure_loss','transition_loss']):
            for h in hs:ax.plot(h.epoch,h[key],color=color,alpha=.15,lw=.7)
            joined=pd.concat(hs).groupby('epoch')[key].mean();ax.plot(joined.index,joined,color=color,label=kind)
    for ax,title in zip(axes,['Total objective','Exposure-rank Huber loss','Auxiliary transition Huber loss']):tidy(ax,title,'Epoch','Training loss')
    axes[0].legend();fig.suptitle('Supplement 1. Optimization diagnostics (not held-out validation)');fig.tight_layout();save(fig,out,'Supplement_1_training')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    rng=np.random.default_rng(0)
    for j,c in enumerate(ORDER):
        f=cohort[cohort.set_membership==c]
        axes[0].scatter(j+rng.uniform(-.2,.2,len(f)),f.normalized_priority,color=COL[c],alpha=.55,s=14)
        axes[1].scatter(j+rng.uniform(-.2,.2,len(f)),f.topk_patient_prevalence,color=COL[c],alpha=.55,s=14)
    for ax,title,y in zip(axes,['Cohort contextual priority','Patient prevalence'],['Normalized median rank priority','Fraction of patients selecting node']):
        ax.set_xticks(range(3),['Set '+c for c in ORDER]);tidy(ax,title,y=y);ax.set_ylim(-.03,1.03)
    fig.suptitle(f'Supplement 2. Node summaries — predetermined dataset {row.seed}');fig.tight_layout();save(fig,out,'Supplement_2_nodes')
    # No claim that a stable node is independently significant or causally exposure-driven.
    chosen=cohort.sort_values(['median_patient_rank','node_idx']).head(min(40,n))
    mat=pn.pivot(index='patient_id',columns='node_idx',values='median_rank').loc[norms.patient_id,chosen.node_idx]
    ordering=np.argsort(norms.exposure,kind='stable');fig,ax=plt.subplots(figsize=(11,5))
    im=ax.imshow((n-mat.to_numpy()[ordering])/(n-1),aspect='auto',vmin=0,vmax=1,cmap='viridis')
    ax.set_xticks(range(len(chosen)),[f'{i} ({c})' for i,c in zip(chosen.node_idx,chosen.set_membership)],rotation=90,fontsize=7)
    tidy(ax,'Supplement 3. Patient-specific contextual priorities','Node index (ground-truth set)','Patients ordered by exposure');fig.colorbar(im,ax=ax,label='Within-patient rank priority')
    fig.tight_layout();save(fig,out,'Supplement_3_patient_nodes')
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    axes[0].scatter(norms.exposure,norms.transformed_L2_squared,s=15)
    tidy(axes[0],'A. Magnitude after loader transformation','Exposure','Total squared paired feature change')
    axes[0].text(.03,.94,f'ρ = {rho(norms.exposure,norms.transformed_L2_squared):.2f}',transform=axes[0].transAxes,va='top')
    sel=audit[audit.measure.isin(['rna_L2_squared','meth_L2_squared','cnv_L2_squared','transformed_L2_squared'])]
    dots(axes[1],sel,'measure','exposure_rho',list(sel.measure.unique()));axes[1].tick_params(axis='x',rotation=35,labelsize=7);axes[1].axhline(0,c='gray',ls='--');axes[1].set_ylim(-1,1)
    tidy(axes[1],'B. Magnitude–exposure audit',y='Spearman correlation per dataset')
    for j,c in enumerate(ORDER):
        d=deg[deg['set']==c].degree;axes[2].plot([j-.15,j+.15],[d.mean()]*2,c=COL[c],lw=6)
    axes[2].set_xticks(range(3),['Set '+c for c in ORDER]);tidy(axes[2],'C. Actual graph connectivity',y='Undirected degree (all nodes equal)')
    fig.suptitle('Supplement 4. Design checks; noisy observed magnitudes need not have exactly zero correlation');fig.tight_layout();save(fig,out,'Supplement_4_design')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    # Average optimization repeats WITHIN each dataset before displaying across-dataset spread.
    rr=runs.groupby(['seed','kind'],as_index=False).A_vs_C_AUROC.mean()
    dots(axes[0],rr,'kind','A_vs_C_AUROC',['initialization','subset','permutation']);axes[0].axhline(.5,c='gray',ls='--');axes[0].set_ylim(0,1)
    tidy(axes[0],'A. Three repetition controls',y='Mean single-run A-versus-C AUROC')
    agg=top.copy();dots(axes[1],agg,'set','recall',ORDER);axes[1].set_ylim(0,1);tidy(axes[1],'B. Recovery adjusted for set size',y='Top-k recall within each set')
    fig.suptitle('Supplement 5. Independent dataset dots; mean ± SD; permutation results are descriptive');fig.tight_layout();save(fig,out,'Supplement_5_robustness')
    (out/'analysis_notes.json').write_text(json.dumps(dict(independent_datasets=len(manifest),representative_dataset=int(row.seed),score='negative median_patient_rank',top_k=a.top_k,uncertainty='Across independent datasets: SD, not node-wise standard error',edge_recovery_evaluated=False,cpg_relationships_evaluated=False,warning='In-sample patient score correlation is an optimization diagnostic, not evidence of recovery. No automatic pass/fail or p-values.'),indent=2))
    print(summary.to_string())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);p.add_argument('--top-k',type=int,default=50)
    main(p.parse_args())
