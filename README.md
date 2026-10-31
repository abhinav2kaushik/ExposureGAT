# ExposureGAT: An Exposure-Aware Graph Attention Framework for Integrative Multi-Omics Data Analysis

## Abstract

ExposureGAT is a graph attention framework for identifying patient-specific, exposure-associated
molecular features and network connections from paired multi-omics data. Each participant is
represented by two graphs on a shared prior network (gene–gene interactions and CpG–gene links): a
reference sample (e.g. adjacent normal tissue) and a paired sample (e.g. tumor), with gene expression,
promoter methylation, copy number and CpG methylation as node features. Both graphs are encoded by the
same attention network; exposure is not a model input and enters only as the training target. From the
within-patient differences, ExposureGAT derives node and edge priorities for every patient and cohort-level
summaries, and quantifies their reproducibility with built-in initialization, patient-subset and
exposure-permutation analyses. Exposure-associated features are identified by calibrating each feature's
alignment with exposure against shuffled-label fits. The framework is evaluated on simulated paired
multi-omics data and on paired TCGA lung adenocarcinoma and squamous cell carcinoma samples.

ExposureGAT is a **descriptive, within-cohort** method. It is not designed to predict exposure in new
individuals. See [Scope and intended use](#scope-and-intended-use).

---

## Quick start

```bash
git clone https://github.com/abhinav2kaushik/ExposureGAT.git
cd ExposureGAT
python3 -m pip install -r requirements.txt
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite"))'

bash exposureGAT.sh demo
```

The demo simulates a small paired dataset (60 patients, 300 genes, 30 CpGs) with planted
exposure-linked genes, fits the model, and writes all tables and figures to `results/demo/`.
It runs in about 2–5 minutes on a laptop CPU.

## Reproducing the TCGA analysis

The processed TCGA-LUAD/LUSC inputs used in the paper (59 paired patients) are available as an archive:

```bash
bash exposureGAT.sh tcga --url <ARCHIVE_URL>
```

This downloads and unpacks the inputs to `tmp/tcga_data/` and runs the full analysis with the
published settings (600 epochs, learning rate 0.003, 10 initializations, 20 patient subsets,
20 exposure permutations). Expect several hours on a CPU. Results are written to `results/tcga/`.

## Running on your own data

```bash
bash exposureGAT.sh run --data path/to/input_folder --name my_study
```

Common options (`bash exposureGAT.sh --help` lists all):

| Option | Default | Meaning |
|---|---|---|
| `--epochs` | 600 | Training epochs (fixed; no early stopping) |
| `--lr` | 0.003 | Learning rate |
| `--hidden` | 8 | Hidden dimension |
| `--initializations` | 10 | Fits with different random seeds |
| `--subsets` | 20 | Fits on random 80% patient subsets |
| `--permutations` | 20 | Fits with shuffled exposure |
| `--exposure_col` | `smoking_intensity` | Exposure column in the metadata |
| `--scaling` | `per_gene` | Feature scaling (`per_gene` or `shared`) |
| `--no-plots` | | Fit only; skip R summaries |

To exclude patients whose exposure is missing or recorded as zero without a confirmed
never-exposed status, prepare the cohort first:

```bash
python3 src/make_analysis_cohort.py --in_dir path/to/inputs --out_dir path/to/analysis_cohort
```

### Input format

One folder containing (rows = samples for feature matrices):

| File | Content |
|---|---|
| `X_rna_gene_features.csv` | samples × genes, expression (e.g. log2-CPM) |
| `X_meth_gene_features.csv` | samples × genes, promoter methylation (beta) |
| `X_cnv_gene_features.csv` | samples × genes, copy number (0 = diploid) |
| `X_cpg_features.csv` | samples × CpGs, beta values |
| `edges_ppi.csv` | `from_idx, to_idx, weight` (gene–gene, both directions) |
| `edges_cpg_to_gene.csv` | `from_idx, to_idx, weight, cpg_id, entrez_id` (CpG → gene) |
| `map_gene_nodes.csv`, `map_cpg_nodes.csv` | identifier → node index |
| `metadata_samples.csv` | `sample_id, patient_id, tissue_type` (`NT`/`TP`), exposure column, optional `smoking_status`, `age`, `gender` |

Gene, methylation and CNV matrices must share gene columns in the same order, and all matrices must
share sample rows. Each patient needs one normal (`NT`) and one tumor (`TP`) sample. Simulated data
may include `sim_ground_truth_nodes.csv` (`node_idx, set_membership`) for recovery metrics.

### Running the model directly

`exposureGAT.sh` calls `exposureGAT.py` and then the R summaries. The model can also be run on its own:

```bash
python3 exposureGAT.py --data_dir path/to/input_folder --out_dir results/my_study/model \
  --epochs 600 --lr 0.003 --initializations 10 --subsets 20 --permutations 20
```

| Argument | Default | Description |
|---|---|---|
| `--data_dir` | required | Input folder (see [Input format](#input-format)) |
| `--out_dir` | required | Output folder; must be new or empty |
| `--exposure_col` | `smoking_intensity` | Exposure column in `metadata_samples.csv` |
| `--scaling` | `per_gene` | Feature scaling: `per_gene` (robust per-gene, fitted on training patients) or `shared` (one scale per modality; used for simulations) |
| `--epochs` | 60 | Training epochs per fit (published TCGA analysis: 600). No early stopping |
| `--lr` | 0.003 | Adam learning rate (weight decay 0.001, gradient clipping 5) |
| `--hidden` | 8 | Hidden dimension of the attention encoder |
| `--initializations` | 10 | Fits on all patients with different random seeds |
| `--subsets` | 20 | Fits on random 80% patient subsets; held-out patients are scored |
| `--permutations` | 20 | Fits with shuffled exposure (descriptive null reference) |
| `--top_k` | 100 | Top nodes per patient used for stability and prevalence summaries |
| `--edge_top_fraction` | 0.10 | Fraction of edges per patient reported as top edges |
| `--batch_size` | 8 | Patients per batch when scoring (training is full-batch) |
| `--seed` | 42 | Base random seed |

Input and output files are described in [FileDescription.md](FileDescription.md).

## Outputs

```
results/<name>/
  model/    raw fitted outputs: one folder per fit, patient/cohort node and edge scores
  csv/      summary tables
    input_audit/          exposure and molecular input checks
    figures/              source data for figures, cohort tables, edge selections
    exposure_associated/  exposure-aligned nodes calibrated against shuffled-label fits
    reproducibility/      stability across seeds and patient subsets
  plots/    figures (PDF and PNG), same subfolders
  logs/     settings and one log per step
```

A full description of every input and output file is in [FileDescription.md](FileDescription.md).

Key tables:

| File | Content |
|---|---|
| `model/patient_node_scores.csv` | per patient and node: priority, rank, stability across seeds/subsets, permutation reference |
| `model/cohort_summary.csv` | cohort node ranking and stability |
| `model/patient_edge_scores.csv.gz`, `model/cohort_edge_scores.csv` | edge priorities (patient and cohort) |
| `csv/exposure_associated/Exposure_associated_node_scores_*.csv` | exposure alignment, null percentile, exposure-associated score per node |
| `csv/reproducibility/Held_out_summary.csv` | fit, held-out agreement and ranking stability |

## Method summary

For patient *i*, normal and tumor graphs are encoded by the same one-layer graph attention encoder
with shared parameters; exposure is **not** a model input. With node embeddings
**h**<sub>*i,v*</sub><sup>N</sup> and **h**<sub>*i,v*</sub><sup>T</sup>, the paired change is
Δ**h**<sub>*i,v*</sub> = **h**<sup>T</sup> − **h**<sup>N</sup>, centered on the cohort mean:
**r**<sub>*i,v*</sub> = Δ**h**<sub>*i,v*</sub> − mean<sub>*j*</sub> Δ**h**<sub>*j,v*</sub>.

- **Node contribution:** *c*<sub>*i,v*</sub> = σ(**w**<sub>g</sub>·**r**<sub>*i,v*</sub>) · (**w**<sub>r</sub>·**r**<sub>*i,v*</sub>) / √*V*
- **Patient exposure score:** *ŝ*<sub>*i*</sub> = *b* + Σ<sub>*v*</sub> *c*<sub>*i,v*</sub>
- **Node priority:** |*c*<sub>*i,v*</sub>|
- **Edge priority:** change in *ŝ*<sub>*i*</sub> when the message along edge *u*→*v* is suppressed in both
  states, with all parameters, other messages and the cohort reference fixed.
- **Training loss:** smooth-L1 between *ŝ*<sub>*i*</sub> and the patient's standardized exposure rank,
  plus 0.05 × reconstruction of paired feature differences, 0.01 × mean node gate and 0.001 × mean
  edge weight.

**Sensitivity analyses** (built in): random initializations, 80% patient subsets (scaling and cohort
reference refitted on training patients; held-out patients scored by the fitted model), and
shuffled-exposure fits used as a descriptive null reference. No node- or edge-level p-values are
computed.

**Exposure-associated nodes** (`src/analysis/exposure_associated_nodes.R`): node priority measures
contribution to the fitted score, not association with exposure. Exposure alignment is the Spearman
correlation between a node's signed contribution and exposure, compared with the pooled alignment
reached in shuffled-label fits; nodes above the 99th percentile are reported as exposure-associated.

## Scope and intended use

- **Descriptive, within-cohort.** Results describe the analyzed cohort. Weak agreement with exposure
  in held-out patients is expected and should not be read as predictive performance.
- **Paired design.** Within-patient differences control patient-level factors shared by both samples
  (age, sex, inherited background). Factors that differ between the paired samples (e.g. histology,
  tumor purity, batch) are not controlled. Exposure effects present equally in both samples (e.g.
  field effects) are outside the method's scope.
- **Exposure values must be real measurements.** Do not code missing exposure as zero; use
  `src/make_analysis_cohort.py`.
- **Priority is not association.** A node's priority reflects its contribution to the fitted score.
  Use the exposure-associated analysis and the shuffled-label comparisons before interpreting
  features as exposure-related.
- **Edges are a functional proxy.** Edge priorities describe how existing connections contribute in
  each patient; they do not establish gained or lost interactions.
- **Covariates.** To account for known covariates, residualize exposure on them and use the residual
  as the exposure column, or analyze strata separately when sample size allows.
- **Training length.** Choose it from the separation between real- and shuffled-label fits before
  examining node or edge results (`src/diagnose_exposure.py --stage opt`), then keep it fixed.

## Repository structure

```
exposureGAT.sh                 single entry point (demo | tcga | run)
exposureGAT.py                 model fitting and sensitivity analyses
src/
  data_loader.py               input loading, pairing and training-only scaling
  simulation_scaling.py        shared scaling used for simulated data
  extended_scores.py           edge scores and patient summaries
  audit_inputs.py              input audit
  make_analysis_cohort.py      exposure-based patient exclusion
  diagnose_exposure.py         training-length and signal diagnostics
  simulation/simulate_data.py  demo / benchmark data simulator
  plotting/plot_results.R      figures and summary tables
  analysis/                    exposure-associated nodes, reproducibility, node heatmaps
  preprocessing/               TCGA input construction, RNA normalization, missingness report
results/  tmp/                 created at run time (not tracked)
```

### Building TCGA inputs from raw GDC data (optional)

`src/preprocessing/build_tcga_inputs.R` builds the inputs from GDC downloads (RNA-seq STAR counts,
450K methylation, ASCAT3 gene-level copy number; the download code is at the top of the script),
`src/preprocessing/normalize_rna.R` applies TMM log2-CPM normalization, and
`src/make_analysis_cohort.py` creates the analysis cohort. These steps require Bioconductor packages
(TCGAbiolinks, SummarizedExperiment, edgeR, limma, biomaRt, STRINGdb, and others listed in each script).
External resources (Ensembl, STRING) are queried at run time, so rebuilding may differ slightly from
the archived inputs; use the archive to reproduce published results exactly.

## Requirements

- Python ≥ 3.9 with numpy, pandas, scipy, torch, torch-geometric (and matplotlib, scikit-learn for diagnostics)
- R ≥ 4.1 with dplyr, tidyr, readr, ggplot2, patchwork, jsonlite
- CPU is sufficient; memory use grows with patients × edges (TCGA run: about 4 GB of output for 50 fits)

## Citation

If you use ExposureGAT, please cite: *[citation to be added]*

## License

*[license to be added]*
