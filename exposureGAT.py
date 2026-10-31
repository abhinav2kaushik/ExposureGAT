"""ExposureGAT: descriptive patient-specific exposure GNN with three repetition controls
(initialization, patient-subset and exposure-permutation fits).
"""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from scipy.stats import rankdata
from data_loader import LungMultiOmicsDataset
from simulation_scaling import SimulationDataset
from extended_scores import export_run, summarize

class Net(nn.Module):
    def __init__(self,n,g,edges,h=8):
        super().__init__();self.g=g
        self.gene=nn.Linear(3,h);self.cpg=nn.Linear(1,h)
        self.selfmap=nn.Linear(h,h,bias=False);self.message=nn.Linear(h,h,bias=False)
        self.key=nn.Linear(2*h,h);self.att=nn.Linear(h,1,bias=False);self.edgegate=nn.Linear(2*h,1)
        self.nodegate=nn.Linear(h,1);self.readout=nn.Linear(h,1,bias=False)
        self.decode=nn.Linear(h,3);self.bias=nn.Parameter(torch.zeros(()))
        self.register_buffer('src',edges[0]);self.register_buffer('dst',edges[1])
        valid=torch.ones(n,3);valid[g:,1:]=0;self.register_buffer('valid',valid)
    def encode(self,x):
        h=torch.tanh(torch.cat([self.gene(x[:,:self.g]),self.cpg(x[:,self.g:,:1])],1))
        pair=torch.cat([h[:,self.src],h[:,self.dst]],-1)
        a=torch.exp(self.att(F.leaky_relu(self.key(pair),.2)).squeeze(-1).clamp(-8,8))
        den=torch.zeros(x.shape[:2]).index_add(1,self.dst,a)
        weights=a/den[:,self.dst].clamp_min(1e-8)*torch.sigmoid(self.edgegate(pair).squeeze(-1))
        agg=torch.zeros_like(h).index_add(1,self.dst,weights[...,None]*self.message(h[:,self.src]))
        base=self.selfmap(h)
        return torch.tanh(base+agg),torch.tanh(base),weights
    def forward(self,n,t,reference=None):
        hn,bn,an=self.encode(n);ht,bt,at=self.encode(t)
        delta=ht-hn
        if reference is None:reference=delta.mean(0)
        r=delta-reference;p=torch.sigmoid(self.nodegate(r)).squeeze(-1)
        c=p*self.readout(r).squeeze(-1)/np.sqrt(r.shape[1])
        local=(bt-bn)-reference
        cl=torch.sigmoid(self.nodegate(local)).squeeze(-1)*self.readout(local).squeeze(-1)/np.sqrt(r.shape[1])
        return self.bias+c.sum(1),c,delta,self.decode(delta)*self.valid,cl,p,(an+at)/2

def rank_target(x):
    r=(rankdata(x,method='average')-.5)/len(x)
    return (r-r.mean())/max(r.std(),1e-6)

def fit(ds,exposure,train,seed,a):
    torch.manual_seed(seed);np.random.seed(seed)
    ds.fit_scalers([ds.patient_ids[i] for i in train])
    normals=[];tumors=[]
    for i in range(len(ds)):
        n,t,_,_=ds[i];normals.append(n.x);tumors.append(t.x)
    n=torch.stack(normals);t=torch.stack(tumors)
    m=Net(ds.n_nodes,ds.n_genes,ds.edge_index,a.hidden)
    y=torch.tensor(rank_target(exposure[train]),dtype=torch.float32)
    opt=torch.optim.Adam(m.parameters(),lr=a.lr,weight_decay=.001)
    scale=(t[train]-n[train]).square().mean((0,1)).sqrt().clamp_min(.03)
    logs=[]
    for epoch in range(a.epochs):
        m.train();pred,c,d,recon,_,gate,edge=m(n[train],t[train])
        le=F.smooth_l1_loss(pred,y)
        err=F.smooth_l1_loss(recon/scale,(t[train]-n[train])/scale,reduction='none')
        lt=(err*m.valid).sum()/(len(train)*m.valid.sum())
        loss=le+.05*lt+.01*gate.mean()+.001*edge.mean()
        opt.zero_grad();loss.backward();nn.utils.clip_grad_norm_(m.parameters(),5);opt.step()
        logs.append(dict(epoch=epoch+1,loss=float(loss.detach()),exposure_loss=float(le.detach()),transition_loss=float(lt.detach())))
    m.eval()
    with torch.no_grad():
        reference=m(n[train],t[train])[2].mean(0)
        scores=[];signed=[];neighbor=[];patient=[]
        for start in range(0,len(ds),a.batch_size):
            pred,c,_,_,local,_,_=m(n[start:start+a.batch_size],t[start:start+a.batch_size],reference)
            signed.append(c.numpy());scores.append(c.abs().numpy());neighbor.append((c-local).numpy());patient.extend(pred.tolist())
    return np.concatenate(scores),np.concatenate(signed),np.concatenate(neighbor),patient,logs,m,reference

def main(a):
    torch.set_num_threads(1)
    if min(a.epochs,a.initializations,a.subsets,a.permutations,a.top_k)<1:raise ValueError('Counts must be positive')
    out=Path(a.out_dir)
    if out.exists() and any(out.iterdir()):raise ValueError('Use a new empty output directory to avoid mixing runs')
    out.mkdir(parents=True,exist_ok=True)
    if not 0<a.edge_top_fraction<=1:raise ValueError('edge_top_fraction must be in (0,1]')
    ds=(SimulationDataset if a.scaling=='shared' else LungMultiOmicsDataset)(a.data_dir,exposure_cols=[a.exposure_col],verbose=False)
    raw=np.array([float(ds.meta[ds.meta.patient_id==p][a.exposure_col].iloc[0]) for p in ds.patient_ids])
    if len(raw)<10 or not np.isfinite(raw).all() or np.ptp(raw)==0:raise ValueError('Need >=10 patients and variable finite exposure')
    allidx=np.arange(len(ds));rng=np.random.default_rng(a.seed);runs={};membership=[];subset_members=[]
    for kind,count in [('initialization',a.initializations),('subset',a.subsets),('permutation',a.permutations)]:
        entries=[]
        for rep in range(count):
            train=np.sort(rng.choice(allidx,max(5,int(.8*len(ds))),replace=False)) if kind=='subset' else allidx
            if kind=='subset':subset_members.append(np.isin(allidx,train))
            labels=rng.permutation(raw) if kind=='permutation' else raw
            # Pair permutation-run seeds with initialization seeds to limit optimization differences.
            seed=a.seed+(rep%a.initializations if kind=='permutation' else rep)
            score,signed,neighbor,patient,history,model,ref=fit(ds,labels,train,seed,a)
            entries.append(score)
            folder=out/f'{kind}_{rep+1}';folder.mkdir(exist_ok=True)
            export_run(model,ds,ref,folder,a.batch_size)
            np.savez_compressed(folder/'attributions.npz',score=score,signed=signed,neighbor_difference=neighbor,patient_score=patient)
            pd.DataFrame(history).to_csv(folder/'history.csv',index=False)
            pd.DataFrame({'patient_id':ds.patient_ids,'included_in_fit':np.isin(allidx,train),'supplied_exposure':labels}).to_csv(folder/'fit_membership.csv',index=False)
            membership.append(dict(kind=kind,repeat=rep+1,seed=seed,n_fit=len(train)))
            if kind=='initialization' and rep==0:
                torch.save(dict(state_dict=model.state_dict(),reference=ref,genes=ds.genes,cpgs=ds.cpgs,config=vars(a),scalers={k:getattr(ds,k) for k in ['_rna_center','_rna_scale','_meth_center','_meth_scale','_cnv_center','_cnv_scale']}),out/'reference_model.pt')
            print(f'{kind} {rep+1}/{count}: completed',flush=True)
        runs[kind]=np.stack(entries)
    k=min(a.top_k,ds.n_nodes)
    def ranks(x):return np.apply_along_axis(lambda v:rankdata(-v,method='average'),-1,x)
    ir=ranks(runs['initialization']);sr=ranks(runs['subset']);pr=ranks(runs['permutation'])
    median_rank=np.median(ir,axis=0)
    # DESCRIPTIVE null comparison, not a calibrated node-wise p-value.
    null_fraction=np.mean(pr<=median_rank[None],axis=0)
    init_freq=(ir<=k).mean(0);sub_freq=(sr<=k).mean(0)
    included=np.stack(subset_members)
    def conditional_frequency(mask):
        numerator=((sr<=k)*mask[:,:,None]).sum(0)
        denominator=mask.sum(0)[:,None]
        return np.divide(numerator,denominator,out=np.full(numerator.shape,np.nan),where=denominator>0)
    freq_in=conditional_frequency(included);freq_out=conditional_frequency(~included)
    names=['GENE:'+g for g in ds.genes]+['CpG:'+c for c in ds.cpgs]
    records=[]
    for i,pid in enumerate(ds.patient_ids):
        for v,label in enumerate(names):records.append(dict(patient_id=pid,node_idx=v,node_label=label,median_priority=float(np.median(runs['initialization'][:,i,v])),median_rank=float(median_rank[i,v]),initialization_topk_frequency=float(init_freq[i,v]),subset_topk_frequency=float(sub_freq[i,v]),subset_included_frequency=float(freq_in[i,v]),subset_excluded_frequency=float(freq_out[i,v]),null_fraction_at_least_as_high_rank=float(null_fraction[i,v])))
    table=pd.DataFrame(records);table.to_csv(out/'patient_node_scores.csv',index=False)
    top=table.sort_values(['patient_id','median_rank','node_idx']).groupby('patient_id',sort=False).head(k);top.to_csv(out/'patient_top_nodes.csv',index=False)
    cohort=table.groupby(['node_idx','node_label']).agg(median_patient_rank=('median_rank','median'),mean_initialization_stability=('initialization_topk_frequency','mean'),mean_subset_stability=('subset_topk_frequency','mean'),mean_null_fraction=('null_fraction_at_least_as_high_rank','mean')).reset_index()
    prevalence=top.groupby('node_idx').patient_id.nunique()/len(ds)
    cohort['topk_patient_prevalence']=cohort.node_idx.map(prevalence).fillna(0);cohort.to_csv(out/'cohort_summary.csv',index=False)
    pd.DataFrame(membership).to_csv(out/'runs.csv',index=False)
    summarize(out,ds,a.edge_top_fraction)
    (out/'config.json').write_text(json.dumps(vars(a),indent=2))
    truth=Path(a.data_dir)/'sim_ground_truth_nodes.csv'
    if truth.exists():
        gt=pd.read_csv(truth)
        merged=cohort.merge(gt[['node_idx','set_membership']],on='node_idx',validate='one_to_one')
        ac=merged[merged.set_membership.isin(['A','C'])]
        positive=(ac.set_membership=='A').to_numpy();values=-ac.median_patient_rank.to_numpy()
        n1=positive.sum();n0=(~positive).sum()
        result={'score':'negative median patient rank','A_vs_C_AUROC':float((rankdata(values)[positive].sum()-n1*(n1+1)/2)/(n1*n0)) if n1*n0 else None}
        (out/'simulation_recovery.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data_dir',required=True);p.add_argument('--out_dir',required=True)
    p.add_argument('--scaling',choices=['shared','per_gene'],default='per_gene');p.add_argument('--exposure_col',default='smoking_intensity')
    for key,val in [('epochs',60),('initializations',10),('subsets',20),('permutations',20),('top_k',100),('hidden',8),('batch_size',8),('seed',42)]:p.add_argument('--'+key,type=int,default=val)
    p.add_argument('--lr',type=float,default=.003)
    p.add_argument('--edge_top_fraction',type=float,default=.10)
    main(p.parse_args())
