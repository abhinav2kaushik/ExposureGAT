"""Read-only audit of EXPORTED TCGA inputs; cannot recover pre-imputation missingness."""
import argparse,json
from pathlib import Path
import pandas as pd
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--data_dir',required=True);p.add_argument('--out_dir',required=True);a=p.parse_args()
d=Path(a.data_dir);o=Path(a.out_dir);o.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(d/'metadata_samples.csv')
cols=[c for c in ['smoking_intensity','smoking_status','age','gender'] if c in m]
conflicts=m.groupby('patient_id')[cols].nunique(dropna=False).gt(1).any(axis=1)
if conflicts.any():raise ValueError('Conflicting paired patient metadata: '+','.join(conflicts[conflicts].index))
pats=m.drop_duplicates('patient_id').copy()
pats['zero_recorded_exposure']=pd.to_numeric(pats.smoking_intensity,errors='coerce').eq(0)
pats['exposure_provenance_status']=np.where(pats.zero_recorded_exposure,'zero requires original-record verification','positive exported value; original provenance still required')
pats.to_csv(o/'exported_patient_audit.csv',index=False)
if 'smoking_status' in pats:pd.crosstab(pats.smoking_status,pats.zero_recorded_exposure).to_csv(o/'smoking_status_by_zero_exposure.csv')
rows=[]
for name in ['rna_gene','meth_gene','cnv_gene','cpg']:
 x=pd.read_csv(d/f'X_{name}_features.csv',index_col=0)
 rows.append(dict(modality=name,samples=x.shape[0],features=x.shape[1],nonfinite_exported_cells=int((~np.isfinite(x.to_numpy(dtype=float))).sum()),
                  constant_sample_rows=int((x.nunique(axis=1,dropna=False)<=1).sum())))
pd.DataFrame(rows).to_csv(o/'exported_molecular_audit.csv',index=False)
x=pd.read_csv(d/'X_cnv_gene_features.csv',index_col=0)
flags=pd.DataFrame({'sample_id':x.index,'all_zero':x.eq(0).all(axis=1).to_numpy(),
                    'constant':x.nunique(axis=1).eq(1).to_numpy()}).merge(m,on='sample_id',validate='one_to_one')
flags[['sample_id','patient_id','tissue_type','all_zero','constant']].to_csv(o/'cnv_profile_flags.csv',index=False)
flags.groupby('tissue_type')[['all_zero','constant']].sum().to_csv(o/'cnv_flags_by_tissue.csv')
summary=dict(patients=len(pats),samples=len(m),zero_exposure=int(pats.zero_recorded_exposure.sum()),
 zero_exposure_status_counts=pats.loc[pats.zero_recorded_exposure,'smoking_status'].value_counts().to_dict() if 'smoking_status' in pats else {},
 unavailable_covariates=[c for c in ['histology','tumor_purity','batch'] if c not in m],
 warning='Counts describe exported inputs, not original clinical or molecular missingness. Do not equate zero with never-smoking. Constant profiles are flags, not proof of imputation.')
(o/'audit_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
