"""Create the analysis cohort: copy an input folder, excluding patients without usable exposure.

A patient is excluded when the exposure value is zero (or missing) and the smoking status
is not a confirmed never-smoker. Confirmed never-smokers are kept with zero exposure.
Input files are never modified; the output folder holds real copies (safe to zip/share).

Example:
  python src/make_analysis_cohort.py --in_dir GNN_Input_Data_v2 --out_dir GNN_Input_Data_v2_clean59
"""
import argparse
import shutil
from pathlib import Path

import pandas as pd

NEVER = 'never|non-smoker|nonsmoker|non smoker|lifelong'


def main(a):
    src, dst = Path(a.in_dir), Path(a.out_dir)
    if dst.exists() and any(dst.iterdir()):
        raise ValueError(f'Output folder must be new or empty: {dst}')
    dst.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(src / 'metadata_samples.csv')
    patients = meta.drop_duplicates('patient_id').copy()
    exposure = pd.to_numeric(patients[a.exposure_col], errors='coerce')
    status = patients[a.status_col].astype(str).str.lower() if a.status_col in patients else pd.Series('', index=patients.index)
    never = status.str.contains(NEVER)
    exclude = (exposure.isna() | (exposure == 0)) & ~never
    excluded = patients.loc[exclude, 'patient_id'].astype(str)
    for f in src.iterdir():
        if f.is_file() and f.name != 'metadata_samples.csv':
            shutil.copy2(f.resolve(), dst / f.name)
    meta[~meta.patient_id.astype(str).isin(excluded)].to_csv(dst / 'metadata_samples.csv', index=False)
    report = patients[['patient_id', a.exposure_col] + ([a.status_col] if a.status_col in patients else [])].copy()
    report['excluded'] = exclude.to_numpy()
    report.to_csv(dst / 'cohort_selection.csv', index=False)
    print(f'Patients: {len(patients)} total, {int(exclude.sum())} excluded '
          f'(zero/missing exposure without never-smoker status), {len(patients) - int(exclude.sum())} retained -> {dst}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--in_dir', required=True)
    p.add_argument('--out_dir', required=True)
    p.add_argument('--exposure_col', default='smoking_intensity')
    p.add_argument('--status_col', default='smoking_status')
    main(p.parse_args())
