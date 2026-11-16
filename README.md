# ExposureGAT: An Exposure-Aware Graph Attention Framework for Integrative Multi-Omics Data Analysis

## Abstract

ExposureGAT is a graph attention framework for identifying patient-specific, exposure-associated
molecular features and network connections from paired multi-omics data. Each participant is
represented by two graphs on a shared prior network (gene–gene interactions and CpG–gene links): a
reference sample (e.g. adjacent normal tissue) and a paired sample (e.g. tumor), with gene expression,
promoter methylation, copy number and CpG methylation as node features. Both graphs are encoded by the
same attention network; exposure is not a model input and enters only as the training target. From the
within-patient differences, ExposureGAT derives a signed contribution of every node in every patient,
node and edge priorities, and the **weighted node score**, its main output: a cohort-level node ranking
by the alignment of each node's patient-specific contributions with exposure. The weighted node score is
produced by every run (simulation, TCGA and your own data) and for every node type (genes, CpGs). Reproducibility and
exposure-label dependence are quantified with built-in initialization, patient-subset and
exposure-permutation analyses. The framework is evaluated on simulated paired multi-omics data and on
paired TCGA lung adenocarcinoma and squamous cell carcinoma samples.

ExposureGAT is a **descriptive, within-cohort** method. It is not designed to predict exposure in new
individuals. See [Scope and intended use](#scope-and-intended-use).

---

## Quick start

```bash
git clone https://github.com/abhinav2kaushik/ExposureGAT.git
cd ExposureGAT
python3 -m pip install -r requirements.txt          # includes Snakemake
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite","ggpubr"))'

snakemake --cores 4 --configfile config/demo.yaml   # quick test, about 2 minutes
```

The demo simulates two small datasets, runs the simulation benchmark and the cohort analysis on
them, and writes every table and figure to `results/demo_simulation/` and `results/demo_cohort/`.
The main result is `results/demo_cohort/weighted_node_score/weighted_node_scores.csv`.

## Workflow

The analysis is a [Snakemake](https://snakemake.readthedocs.io) workflow (`Snakefile`); all settings
are in `config/config.yaml`.

| Command | What it does |
|---|---|
| `snakemake --cores 8` | everything in the manuscript: TCGA analysis and simulation benchmark |
| `snakemake --cores 8 cohort` | cohort analysis only (TCGA; see [TCGA.md](TCGA.md)) |
| `snakemake --cores 8 simulation` | simulation benchmark only |
| `snakemake --cores 4 --configfile config/demo.yaml` | quick test on small simulated data |
| `snakemake -n` | dry run: list the steps that would run |

Snakemake runs only the steps whose outputs are missing or out of date, and runs independent steps in
parallel (for example, the ten simulated datasets). Each step writes a log to `results/<name>/logs/`.

**Cohort analysis** (`results/<cohort name>/`):

| Step (rule) | Output |
|---|---|
| `audit_inputs` | `input_audit/`: input checks |
| `fit_model` | `model/`: 10 initializations, 20 fits on random 80% patient subsets, 20 fits with permuted exposure |
| `weighted_node_score` | `weighted_node_score/` (**main output**); Figure 3B–C |
| `label_controls` | `label_controls/`: real- vs shuffled-label statistics; Supplementary Figures 10–11 |
| `held_out_summary` | `reproducibility/`: held-out agreement and ranking stability |
| `direction_of_change`, `marker_ranks` | `genes/`: measured direction of change of the top genes; marker positions |
| `gene_signatures` | `gene_signatures/`: MSigDB signature representation; Supplementary Figures 9, 12, 13 |
| `cohort_figures` | Figure 3A; Supplementary Figures 6–8 |

**Simulation benchmark** (`results/<simulation name>/`):

| Step (rule) | Output |
|---|---|
| `simulate_dataset` | `data/seed_<n>/`: one independent simulated dataset per seed |
| `fit_simulated`, `weighted_simulated` | `seed_<n>/model/`, `seed_<n>/weighted_node_score/` |
| `simulation_tables` | `tables/`: AUROCs, top-k recovery, repetition controls, design checks |
| `simulation_figures` | Figure 2B–D; Supplementary Figures 1–5 |

All figures are written to `figures/` as PDF and PNG, named after the manuscript figure they
correspond to (e.g. `Figure_3A.pdf`, `Supplementary_Figure_07.pdf`); figure data are saved next to
them as CSV. Figure 1 and Figure 2A are schematics and are not generated.

## Simulation design

The simulator (`src/simulation/simulate_three_set.py`) creates paired data on a regular graph (every
node has the same degree) with three planted gene sets: **A**, whose paired changes are redistributed
according to exposure while the total planted change is the same in every patient; **C**, with the same
changes assigned independently of exposure; and **B**, background. Distinguishing A from C tests whether
exposure-associated change is ranked above generic paired change. Settings (seeds, patients, genes, set
sizes, degree, fits per dataset) are in the `simulation` section of `config/config.yaml`; the published
benchmark is 10 datasets × 50 fits (several hours on a CPU).

## Running on TCGA data

Follow **[TCGA.md](TCGA.md)**: download the processed inputs (Harvard Dataverse,
https://doi.org/10.7910/DVN/ABBWZG) into the repository and run `snakemake --cores 8 cohort`.

## Running on your own data

Copy the configuration, point the `cohort` section to your input folder and run the cohort workflow:

```bash
cp config/config.yaml config/my_study.yaml        # edit cohort: name, data_dir, exposure_col, ...
snakemake --cores 8 cohort --configfile config/my_study.yaml
```

| Setting (`cohort:`) | Published value | Meaning |
|---|---|---|
| `name` | `tcga` | results go to `results/<name>/` |
| `data_dir` | `data/tcga` | input folder (format below) |
| `exposure_col` | `smoking_intensity` | exposure column in `metadata_samples.csv` |
| `reference_label`, `transition_label` | `NT`, `TP` | `tissue_type` values of the two paired states |
| `epochs`, `lr`, `hidden` | 600, 0.003, 8 | training (fixed length; no early stopping) |
| `initializations`, `subsets`, `permutations` | 10, 20, 20 | repetition analyses |
| `top_k` | 100 | top nodes per patient for stability summaries |
| `gamma`, `top_fraction` | 0.10, 0.10 | weighted node score settings |
| `gmt_dir`, `catalog` | `MSigDB`, `resources/signature_catalog.csv` | gene-signature step (set `gmt_dir: ""` to skip) |
| `marker_genes`, `marker_cpgs` | TCGA smoking anchors | Entrez IDs / CpG probes to locate in the rankings (`""` to skip) |
| `gene_symbols` | `""` | optional CSV `entrez_id,symbol` for gene labels in Figure 3B |
| `scaling`, `seed` | `per_gene`, 42 | feature scaling; base random seed |

Exclude patients whose exposure is missing (do not code it as zero) before running:

```bash
python3 src/make_analysis_cohort.py --in_dir path/to/inputs --out_dir path/to/analysis_cohort
```

### Input format

One folder containing (rows = samples for feature matrices):

| File | Content |
|---|---|
| `X_rna_gene_features.csv` | samples × genes, expression (e.g. TMM log2-CPM) |
| `X_meth_gene_features.csv` | samples × genes, promoter methylation (beta) |
| `X_cnv_gene_features.csv` | samples × genes, copy number (0 = diploid) |
| `X_cpg_features.csv` | samples × CpGs, beta values |
| `edges_ppi.csv` | `from_idx, to_idx, weight` (gene–gene, both directions) |
| `edges_cpg_to_gene.csv` | `from_idx, to_idx, weight, cpg_id, entrez_id` (CpG → gene) |
| `map_gene_nodes.csv`, `map_cpg_nodes.csv` | identifier → node index |
| `metadata_samples.csv` | `sample_id, patient_id, tissue_type` (`NT`/`TP`), exposure column, optional `smoking_status`, `age`, `gender` |

Gene, methylation and CNV matrices share gene columns in the same order, and all matrices share
sample rows. Each patient needs one reference (`NT`) and one transition (`TP`) sample. Simulated data
include `sim_ground_truth_nodes.csv` (`node_idx, set_membership`) for recovery metrics. Every input and
output file is described in [FileDescription.md](FileDescription.md).

## Main output

| File | Content |
|---|---|
| `weighted_node_score/weighted_node_scores.csv` | **main output**: genes ranked by weighted node score, with exposure alignment and paired-transition magnitude |
| `weighted_node_score/weighted_node_scores_CpG.csv` | the same for CpG nodes, ranked among CpGs (when the graph has CpG nodes) |
| `weighted_node_score/patient_weighted_contributions.csv.gz` | top genes per patient by weighted contribution (`_CpG` file for CpGs) |
| `label_controls/statistical_comparisons.csv` | real- versus shuffled-label comparisons (run-group randomization, Holm-adjusted) |
| `reproducibility/Held_out_summary.csv` | in-sample fit, held-out agreement and ranking stability |
| `model/patient_node_scores.csv`, `model/patient_edge_scores.csv.gz` | node and edge priorities per patient |

Every output file is described in [FileDescription.md](FileDescription.md).

## Method summary

For patient *p*, the reference and transition graphs are encoded by the same one-layer graph attention
encoder; exposure is **not** a model input. The paired embedding change of node *i* is
Δ**h**(*p,i*) = **h**<sup>T</sup>(*p,i*) − **h**<sup>N</sup>(*p,i*), centered on the cohort mean:
**r**(*p,i*) = Δ**h**(*p,i*) − mean<sub>*q*</sub> Δ**h**(*q,i*).

- **Signed contribution:** *c*(*p,i*) = σ(**w**<sub>n</sub>·**r** + *b*<sub>n</sub>) · (**w**<sub>o</sub>·**r**) / √*V*
- **Patient exposure score:** *ŝ*(*p*) = *b* + Σ<sub>*i*</sub> *c*(*p,i*), trained (smooth-L1) to fit the
  patient's standardized exposure rank, plus 0.05 × reconstruction of the paired feature differences,
  0.01 × mean node gate and 0.001 × mean gated attention.
- **Node priority:** |*c*(*p,i*)|, ranked within each patient.
- **Edge priority:** absolute change in *ŝ*(*p*) when the message along an edge is suppressed in both
  states, with all parameters, other messages and the cohort reference fixed.
- **Weighted node score** (`src/scoring/weighted_node_score.py`): for each node, exposure alignment
  *A* is the Spearman correlation across patients between *c*(*p,i*) and exposure; *M* is the median
  paired-transition magnitude ‖Δ**h**‖ and *T* its percentile among nodes of the same type (genes, CpGs);
  *W* = max(*A*, 0) × (0.9 + 0.1 *T*). The cohort score uses the medians of *A* and *M* over the
  real-label initializations. The patient-level weighted contribution is *W* × |*c*(*p,i*)|.

Node priority reflects a node's contribution to the fitted score, which can be driven by large paired
change regardless of exposure; exposure-related interpretation therefore uses the weighted node score.
No node- or edge-level P values are computed; P values are reported only for aggregate comparisons
between groups of fits.

## Scope and intended use

- **Descriptive, within-cohort.** Results describe the analyzed cohort. Weak agreement with exposure
  in held-out patients is expected and should not be read as predictive performance.
- **Paired design.** Pairing reduces baseline variation shared by both samples (e.g. related to age,
  sex, inherited background) but does not establish independence from factors that modify the paired
  difference or differ between the paired samples (e.g. histology, tumor purity, batch). Effects shared
  equally by both samples (e.g. field effects) are outside the method's scope.
- **Exposure values must be real measurements.** Do not code missing exposure as zero; use
  `src/make_analysis_cohort.py`.
- **The weighted node score uses the exposure labels.** Its relation to exposure is therefore not
  independent validation; compare with the shuffled-label fits (`label_controls/`).
- **Edges are a functional proxy.** Edge priorities describe how existing connections contribute in each
  patient; they do not establish changes in biological interactions.
- **Covariates.** To account for known covariates, residualize exposure on them and use the residual
  as the exposure column, or analyze strata separately when sample size allows.
- **Training length.** Choose it from the separation between real- and shuffled-label fits before
  examining node or edge results (`src/diagnose_exposure.py --stage opt`), then keep it fixed.

## Repository structure

```
Snakefile                            the workflow (run with snakemake)
config/config.yaml                   published settings (cohort and simulation)
config/demo.yaml                     small settings for a quick test
resources/signature_catalog.csv      smoking-, nicotine- and lung cancer-related MSigDB signatures (49 retained)
MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip   MSigDB v2026.1 gene sets (Entrez IDs), read directly by the gene-signature step
src/
  exposuregat_model.py               model fitting, scoring and repetition analyses
  data_loader.py, simulation_scaling.py   input loading, pairing and scaling
  extended_scores.py                 edge priorities and patient summaries
  scoring/weighted_node_score.py     weighted node score (main output); Figure 3B-C
  analysis/                          label controls, gene signatures, direction of change, marker ranks, held-out summary
  figures/                           cohort_figures.R (Figure 3A, Supp. Figs 6-8), label_control_figures.R (10-11),
                                     gene_signature_figures.R (9, 12, 13), simulation_figures.R (Figure 2B-D, 1-5)
  simulation/                        simulate_three_set.py, simulation_summary.py
  audit_inputs.py, make_analysis_cohort.py, diagnose_exposure.py
  preprocessing/                     TCGA input construction, RNA normalization, missingness report
data/  results/                      created at run time (not tracked)
```

### Building TCGA inputs from raw GDC data (optional)

`src/preprocessing/build_tcga_inputs.R` builds the inputs from GDC downloads (RNA-seq STAR counts,
450K methylation, ASCAT3 gene-level copy number; the download code is at the top of the script),
`src/preprocessing/normalize_rna.R` applies TMM log2-CPM normalization, and
`src/make_analysis_cohort.py` creates the analysis cohort. These steps require Bioconductor packages
(TCGAbiolinks, SummarizedExperiment, edgeR, limma, biomaRt, STRINGdb and others listed in each script).
External resources (Ensembl, STRING) are queried at run time, so rebuilding may differ slightly from
the archived inputs; use the archive to reproduce published results exactly.

## Requirements

- Python ≥ 3.9 with snakemake, numpy, pandas, scipy, torch, torch-geometric, matplotlib (and scikit-learn for diagnostics)
- R ≥ 4.1 with dplyr, tidyr, readr, ggplot2, patchwork, jsonlite, ggpubr
- CPU is sufficient; memory use grows with patients × edges (TCGA run: about 4 GB of output for 50 fits)

Versions used for the published analysis: Python 3.11.11, PyTorch 2.3.0, PyTorch Geometric 2.7.0,
NumPy 1.26.4, pandas 2.3.3, SciPy 1.16.3; R 4.4.0. The workflow was tested with Snakemake 9. If `pip install torch-geometric` fails, install
PyTorch first and follow https://pytorch-geometric.readthedocs.io for your platform.

**Reproducibility.** All fits are seeded (`seed` in the configuration, default 42), but floating-point results can differ
slightly across operating systems, CPU types and library versions. Reruns therefore reproduce the
published values closely but not always to the last digit; the published outputs were produced with
the versions above.

**Gene sets.** `MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip` is MSigDB v2026.1 (https://www.gsea-msigdb.org),
distributed under the MSigDB license terms; some collections (for example KEGG) carry their own terms
of use.

## Citation

If you use ExposureGAT, please cite: Kaushik A, Johnson M. ExposureGAT: An Exposure-Aware Graph Attention
Framework for Integrative Multi-Omics Data Analysis. *[journal and year to be added]*

## License

*[license to be added]*
