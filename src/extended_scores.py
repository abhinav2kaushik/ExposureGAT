"""Post-fit diagnostics for the unchanged one-layer personalized GNN.
Exact local message suppression, not retraining or attention renormalization.
"""
import gzip,json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata

FIELDS=['signed_exposure_contribution','exposure_priority','message_transition_norm',
        'reference_attention','transition_attention','reference_gated_attention',
        'transition_gated_attention','attention_change','gated_attention_change']

def contextual_contribution(model,r):
    return (torch.sigmoid(model.nodegate(r)).squeeze(-1)*
            model.readout(r).squeeze(-1)/np.sqrt(model.valid.shape[0]))

@torch.no_grad()
def inspect_batch(model,normal,tumor,reference,edge_chunk=2048):
    """For one message-passing layer, edge suppression changes only its destination.
    Both state messages are removed together, with all fitted parameters, other
    messages, attention normalization and cohort reference held fixed.
    """
    states=[]
    for x in (normal,tumor):
        h=torch.tanh(torch.cat([model.gene(x[:,:model.g]),model.cpg(x[:,model.g:,:1])],1))
        pair=torch.cat([h[:,model.src],h[:,model.dst]],-1)
        a=torch.exp(model.att(torch.nn.functional.leaky_relu(model.key(pair),.2)).squeeze(-1).clamp(-8,8))
        den=torch.zeros(x.shape[:2],device=x.device).index_add(1,model.dst,a)
        alpha=a/den[:,model.dst].clamp_min(1e-8)
        weight=alpha*torch.sigmoid(model.edgegate(pair).squeeze(-1))
        msg=weight[...,None]*model.message(h[:,model.src])
        pre=model.selfmap(h)+torch.zeros_like(h).index_add(1,model.dst,msg)
        states.append((pre,msg,alpha,weight))
    pn,mn,an,wn=states[0];pt,mt,at,wt=states[1]
    d=torch.tanh(pt)-torch.tanh(pn)
    c=contextual_contribution(model,d-reference)
    signed=torch.empty(normal.shape[0],len(model.src),device=normal.device)
    for start in range(0,len(model.src),edge_chunk):
        sl=slice(start,start+edge_chunk);dst=model.dst[sl]
        suppressed=torch.tanh(pt[:,dst]-mt[:,sl])-torch.tanh(pn[:,dst]-mn[:,sl])
        without=contextual_contribution(model,suppressed-reference[dst])
        signed[:,sl]=c[:,dst]-without
    edge=dict(signed_exposure_contribution=signed,exposure_priority=signed.abs(),
              message_transition_norm=(mt-mn).norm(dim=-1),reference_attention=an,
              transition_attention=at,reference_gated_attention=wn,
              transition_gated_attention=wt,attention_change=at-an,
              gated_attention_change=wt-wn)
    node_transition=d.norm(dim=-1)
    patient=dict(predicted_exposure_rank_score=model.bias+c.sum(1),
                 transition_score=node_transition.mean(1),
                 exposure_contribution_magnitude=c.abs().sum(1),
                 scaled_input_transition_rms=((tumor-normal).square()*model.valid).sum((1,2)).div(model.valid.sum()).sqrt())
    return {k:v.cpu().numpy() for k,v in edge.items()},node_transition.cpu().numpy(),{k:v.cpu().numpy() for k,v in patient.items()}


def export_run(model,ds,reference,folder,batch_size):
    folder=Path(folder)
    # Memory-mapped structured array: downstream summaries read one patient at a
    # time instead of stacking all repetitions into RAM.
    target=np.lib.format.open_memmap(folder/'edge_attributions.npy',mode='w+',
              dtype=[(k,'<f4') for k in FIELDS],shape=(len(ds),len(model.src)))
    node=np.lib.format.open_memmap(folder/'node_transition.npy',mode='w+',dtype='float32',shape=(len(ds),ds.n_nodes))
    rows=[]
    for start in range(0,len(ds),batch_size):
        samples=[ds[i] for i in range(start,min(start+batch_size,len(ds)))]
        normal=torch.stack([x[0].x for x in samples]);tumor=torch.stack([x[1].x for x in samples])
        edges,transition,patients=inspect_batch(model,normal,tumor,reference)
        stop=start+len(samples)
        for name,values in edges.items():target[name][start:stop]=values
        node[start:stop]=transition
        rows.append(pd.DataFrame(dict(patient_id=ds.patient_ids[start:stop],**patients)))
    target.flush();node.flush()
    pd.concat(rows,ignore_index=True).to_csv(folder/'patient_scores.csv',index=False)


def summarize(out,ds,edge_top_fraction=.10):
    out=Path(out);runs=pd.read_csv(out/'runs.csv')
    labels=['GENE:'+g for g in ds.genes]+['CpG:'+c for c in ds.cpgs]
    src,dst=ds.edge_index.numpy();e=len(src);n=len(ds);v=ds.n_nodes
    if not e:raise ValueError('Cannot summarize an empty graph')
    edge_ids=np.arange(e);k=max(1,int(np.ceil(edge_top_fraction*e)))
    index=pd.DataFrame(dict(edge_id=edge_ids,source_idx=src,target_idx=dst,
          source_label=[labels[x] for x in src],target_label=[labels[x] for x in dst],
          edge_type=['gene_to_gene' if s<ds.n_genes and t<ds.n_genes else 'cpg_to_gene' if s>=ds.n_genes and t<ds.n_genes else 'gene_to_cpg' if s<ds.n_genes else 'cpg_to_cpg' for s,t in zip(src,dst)]))
    index.to_csv(out/'edge_index.csv',index=False)
    arrays={};members={};patient_runs=[];nodearrays=[]
    for r in runs.itertuples():
        folder=out/f'{r.kind}_{r.repeat}'
        arrays.setdefault(r.kind,[]).append(np.load(folder/'edge_attributions.npy',mmap_mode='r'))
        m=pd.read_csv(folder/'fit_membership.csv');assert list(m.patient_id)==list(ds.patient_ids)
        members.setdefault(r.kind,[]).append(m.included_in_fit.to_numpy(bool))
        p=pd.read_csv(folder/'patient_scores.csv');p['kind']=r.kind;p['repeat']=r.repeat
        p['included_in_fit']=m.included_in_fit;p['supplied_exposure']=m.supplied_exposure
        patient_runs.append(p)
        if r.kind=='initialization':nodearrays.append(np.load(folder/'node_transition.npy',mmap_mode='r'))
    allpatients=pd.concat(patient_runs,ignore_index=True);allpatients.to_csv(out/'patient_scores_by_run.csv',index=False)
    measures=['predicted_exposure_rank_score','transition_score','exposure_contribution_magnitude','scaled_input_transition_rms']
    init=allpatients[allpatients.kind=='initialization']
    summary=init.groupby('patient_id')[measures].agg(['median','min','max','std'])
    summary.columns=['_'.join(c) for c in summary.columns];summary.reset_index().to_csv(out/'patient_scores.csv',index=False)
    aggregate={key:np.empty((n,e),dtype=np.float32) for key in ['median_rank','median_priority','median_message_transition_norm','initialization_topk_frequency','subset_topk_frequency','null_fraction_at_least_as_high_rank']}
    selected=np.zeros(e);node_transition=np.empty((n,v),dtype=np.float32)
    degree_rows=[];toprows=[]
    # Large patient-by-edge table is compressed and streamed.
    with gzip.open(out/'patient_edge_scores.csv.gz','wt') as fh:
        for i,pid in enumerate(ds.patient_ids):
            ranks={kind:np.stack([rankdata(-a['exposure_priority'][i],method='average') for a in ar]) for kind,ar in arrays.items()}
            medrank=np.median(ranks['initialization'],axis=0)
            metrics={f'median_{key}':np.median(np.stack([a[key][i] for a in arrays['initialization']]),axis=0) for key in FIELDS}
            metrics['median_priority']=metrics.pop('median_exposure_priority')
            metrics.update(median_rank=medrank,initialization_topk_frequency=(ranks['initialization']<=k).mean(0),
                subset_topk_frequency=(ranks['subset']<=k).mean(0),
                null_fraction_at_least_as_high_rank=(ranks['permutation']<=medrank).mean(0))
            for included,name in [(True,'subset_included_frequency'),(False,'subset_excluded_frequency')]:
                mask=np.array([x[i]==included for x in members['subset']])
                metrics[name]=(ranks['subset'][mask]<=k).mean(0) if mask.any() else np.full(e,np.nan)
            for key in aggregate:aggregate[key][i]=metrics[key]
            frame=index.copy();frame.insert(0,'patient_id',pid)
            for key,val in metrics.items():frame[key]=val
            frame.to_csv(fh,index=False,header=i==0)
            top=np.lexsort((edge_ids,medrank))[:k];selected[top]+=1
            toprows.append(frame.iloc[top])
            # Degree counts directed message channels; incoming/outgoing stay separate.
            incoming=np.bincount(dst[top],minlength=v);outgoing=np.bincount(src[top],minlength=v)
            instrength=np.bincount(dst,weights=metrics['median_priority'],minlength=v)
            outstrength=np.bincount(src,weights=metrics['median_priority'],minlength=v)
            degree_rows.append(pd.DataFrame(dict(patient_id=pid,node_idx=np.arange(v),node_label=labels,
                selected_incoming_degree=incoming,selected_outgoing_degree=outgoing,
                incoming_exposure_strength=instrength,outgoing_exposure_strength=outstrength)))
            node_transition[i]=np.median(np.stack([a[i] for a in nodearrays]),axis=0)
    pd.concat(toprows,ignore_index=True).to_csv(out/'patient_top_edges.csv.gz',index=False,compression='gzip')
    cohort=index.copy()
    cohort['median_patient_rank']=np.median(aggregate['median_rank'],axis=0)
    cohort['cohort_exposure_priority']=-cohort.median_patient_rank
    cohort['median_patient_priority']=np.median(aggregate['median_priority'],axis=0)
    cohort['median_patient_message_transition_norm']=np.median(aggregate['median_message_transition_norm'],axis=0)
    for col,new in [('initialization_topk_frequency','mean_initialization_stability'),('subset_topk_frequency','mean_subset_stability'),('null_fraction_at_least_as_high_rank','mean_null_fraction')]:cohort[new]=aggregate[col].mean(0)
    cohort['topk_patient_prevalence']=selected/n
    cohort.to_csv(out/'cohort_edge_scores.csv',index=False)
    degree=pd.concat(degree_rows,ignore_index=True);degree.to_csv(out/'patient_node_edge_summary.csv.gz',index=False,compression='gzip')
    degree.groupby(['node_idx','node_label'],as_index=False)[['selected_incoming_degree','selected_outgoing_degree','incoming_exposure_strength','outgoing_exposure_strength']].mean().to_csv(out/'cohort_node_edge_summary.csv',index=False)
    pd.DataFrame(dict(patient_id=np.repeat(ds.patient_ids,v),node_idx=np.tile(np.arange(v),n),node_label=labels*n,transition_score=node_transition.ravel())).to_csv(out/'patient_node_transition_scores.csv.gz',index=False,compression='gzip')
    pd.DataFrame(dict(node_idx=np.arange(v),node_label=labels,median_patient_transition_score=np.median(node_transition,axis=0),mean_patient_transition_score=node_transition.mean(0))).to_csv(out/'cohort_node_transition_scores.csv',index=False)
    (out/'extended_score_metadata.json').write_text(json.dumps(dict(version=1,edge_count=e,edge_top_fraction=edge_top_fraction,edge_top_k=k,
         edge_orientation='directed model message channel; reciprocal edges remain separate',
         suppression='both states; fixed attention normalization, other messages, parameters and fit-cohort reference',
         edge_priority='absolute fitted exposure-score difference upon suppressing this edge message',
         transition_score='mean node paired contextual embedding L2 change, without cohort centering',
         null_comparison='descriptive single-permutation rank versus observed initialization-consensus rank; not a p-value',
         edge_ground_truth_validated=False),indent=2))
