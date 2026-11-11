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
produced by every run (demo, simulation, TCGA and your own data) and for every node type (genes, CpGs). Reproducibility and
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
python3 -m pip install -r requirements.txt
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite","ggpubr"))'

python3 exposureGAT.py demo
```

The demo simulates a small paired dataset (60 patients, 80 genes) with the manuscript's simulation
design, fits the model with a few repetitions, and writes all tables and figures to `results/demo/`.
It takes about 1–2 minutes on a laptop CPU. The main result is
`results/demo/csv/weighted_node_score/weighted_node_scores.csv`.

## Commands

| Command | What it does |
|---|---|
| `python3 exposureGAT.py demo` | small simulated dataset, full pipeline |
| `python3 exposureGAT.py simulation` | the simulation benchmark of the manuscript (10 datasets; Figure 2, Supplementary Figures 1–5) |
| `python3 exposureGAT.py tcga` | the published TCGA analysis (see [TCGA.md](TCGA.md)) |
| `python3 exposureGAT.py run --data DIR` | your own data |

Every analysis runs the same steps and writes them to `results/<name>/`:

1. **Input audit**
2. **Model fitting**: 10 random initializations on all patients, 20 fits on random 80% patient
   subsets and 20 fits with permuted exposure (published settings: 600 epochs, learning rate 0.003,
   hidden dimension 8).
3. **Weighted node score** (main output): consensus ranking of gene nodes (and, separately, of CpG
   nodes), per-fit scores, patient-level weighted contributions, top-node heatmaps and
   exposure-alignment summaries.
4. **Exposure-label controls**: agreement between real-label fits and between real- and shuffled-label
   fits, number of exposure-aligned genes, and fixed top-5% agreement.
5. **Direction of molecular change** of the top genes, and positions of marker genes (`--markers`).
6. **Gene-signature representation** (optional, with `--gmt`): predefined MSigDB signatures among the
   top genes under real and shuffled labels.
7. **Model figures and reproducibility** (fit diagnostics, node and edge priorities, held-out and
   ranking stability).

## Running on simulated data

The simulator (`src/simulation/simulate_three_set.py`) creates paired data on a regular graph (every
node has the same degree) with three planted gene sets: **A**, whose paired changes are redistributed
according to exposure while the total planted change is the same in every patient; **C**, with the same
changes assigned independently of exposure; and **B**, background. Distinguishing A from C tests whether
exposure-associated change is ranked above generic paired change.

Reproduce the manuscript benchmark (500 model fits; several hours on a CPU):

```bash
python3 exposureGAT.py simulation
```

| Option | Default | Meaning |
|---|---|---|
| `--seeds` | 1–10 | independent datasets (graph, set assignment, exposure and noise resampled) |
| `--patients` / `--genes` | 250 / 400 | subjects and gene nodes per dataset |
| `--set_a` / `--set_c` | 80 / 80 | genes in sets A and C (the rest are set B) |
| `--degree` | 8 | degree of every node |
| `--epochs`, `--initializations`, `--subsets`, `--permutations` | 60, 10, 20, 20 | fits per dataset |

Results are written to `results/simulation/`: one model and one weighted node score folder per dataset
(`seed_<n>/model/`, `seed_<n>/weighted_node_score/`), `figures/` (Figure 2,
`AUROC_summary.csv`, `AUROC_by_dataset.csv`, `repetition_AUROC.csv`, `topk_recovery.csv`) and
`supplementary_figures/` (Supplementary Figures 1–5). A quick check with small settings:

```bash
python3 exposureGAT.py simulation --name sim_check --seeds 1 2 --patients 60 --genes 80 \
  --set_a 16 --set_c 16 --degree 6 --epochs 20 --initializations 2 --subsets 2 --permutations 2
```

## Running on TCGA data

Follow **[TCGA.md](TCGA.md)**: download the processed inputs (Harvard Dataverse,
https://doi.org/10.7910/DVN/ABBWZG), run `python3 exposureGAT.py tcga`, and regenerate tables and figures.

## Running on your own data

```bash
python3 exposureGAT.py run --data path/to/input_folder --name my_study
```

| Option | Default | Meaning |
|---|---|---|
| `--epochs` | 600 | training epochs (fixed; no early stopping) |
| `--lr` | 0.003 | learning rate |
| `--hidden` | 8 | hidden dimension |
| `--initializations` | 10 | fits with different random seeds |
| `--subsets` | 20 | fits on random 80% patient subsets |
| `--permutations` | 20 | fits with shuffled exposure |
| `--exposure_col` | `smoking_intensity` | exposure column in the metadata |
| `--top_k` | 100 | top nodes per patient used for stability summaries |
| `--gamma` | 0.10 | weight of the paired-transition percentile in the weighted node score |
| `--top_fraction` | 0.10 | fraction of genes reported as top genes |
| `--gmt` | | folder with MSigDB GMT files (Entrez IDs, `.gmt` or zipped `.gmt.zip`; the repository includes `MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip`) for the gene-signature step |
| `--catalog` | `resources/signature_catalog.csv` | signatures used by the gene-signature step |
| `--markers`, `--marker_cpgs` | | Entrez IDs / CpG probes to locate in the rankings |
| `--reference_label`, `--transition_label` | `NT`, `TP` | `tissue_type` values of the two paired states |
| `--scaling` | `per_gene` | `per_gene` (robust per gene, fitted on the fitted patients) or `shared` (simulations) |
| `--seed` | 42 | base random seed |
| `--name` | data folder name | results go to `results/NAME` |
| `--no-plots` | | write tables only |
| `--rscript` | `Rscript` | path to Rscript, if not on the PATH |

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

## Outputs

```
results/<name>/
  model/    fitted outputs: one folder per fit, patient and cohort node and edge priorities
  csv/      tables          plots/  figures (PDF and PNG)          logs/  settings and step logs
    weighted_node_score/  main output: weighted node scores, patient weighted contributions
    label_controls/       real- versus shuffled-label agreement, aligned-gene counts, statistics
    genes/                direction of molecular change of top genes, marker positions
    gene_signatures/      MSigDB signature representation (with --gmt)
    figures/              model diagnostics, node and edge priorities
    reproducibility/      held-out agreement and ranking stability
    input_audit/          input checks
```

Key tables:

| File | Content |
|---|---|
| `csv/weighted_node_score/weighted_node_scores.csv` | **main output**: genes ranked by weighted node score, with exposure alignment and paired-transition magnitude |
| `csv/weighted_node_score/weighted_node_scores_CpG.csv` | the same for CpG nodes, ranked among CpGs (when the graph has CpG nodes) |
| `csv/weighted_node_score/patient_weighted_contributions.csv.gz` | top genes per patient by weighted contribution (`_CpG` file for CpGs) |
| `csv/label_controls/statistical_comparisons.csv` | real- versus shuffled-label comparisons (run-group randomization, Holm-adjusted) |
| `csv/genes/direction_of_change_top_genes.csv` | measured RNA and methylation change of the top genes |
| `csv/gene_signatures/gene_signature_table.csv` | signature representation under real and shuffled labels |
| `csv/reproducibility/Held_out_summary.csv` | in-sample fit, held-out agreement and ranking stability |
| `model/patient_node_scores.csv`, `model/patient_edge_scores.csv.gz` | node and edge priorities per patient |

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
  independent validation; compare with the shuffled-label fits (`csv/label_controls/`).
- **Edges are a functional proxy.** Edge priorities describe how existing connections contribute in each
  patient; they do not establish changes in biological interactions.
- **Covariates.** To account for known covariates, residualize exposure on them and use the residual
  as the exposure column, or analyze strata separately when sample size allows.
- **Training length.** Choose it from the separation between real- and shuffled-label fits before
  examining node or edge results (`src/diagnose_exposure.py --stage opt`), then keep it fixed.

## Repository structure

```
exposureGAT.py                       the only command to run: demo | simulation | tcga | run
resources/signature_catalog.csv      smoking-, nicotine- and lung cancer-related MSigDB signatures (49 retained)
MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip   MSigDB v2026.1 gene sets (Entrez IDs), read directly by the gene-signature step
src/
  exposuregat_model.py               model fitting, scoring and repetition analyses
  data_loader.py, simulation_scaling.py   input loading, pairing and scaling
  extended_scores.py                 edge priorities and patient summaries
  scoring/weighted_node_score.py     weighted node score (main output)
  analysis/label_controls.py         real- versus shuffled-label controls and statistics
  analysis/gene_signatures.py        MSigDB signature representation
  analysis/direction_of_change.py    measured direction of change of top genes
  analysis/marker_ranks.py           positions of marker genes and CpGs
  analysis/held_out_summary.R        held-out agreement and ranking stability
  plotting/                          plot_results.R, plot_label_controls.R, plot_gene_signatures.R
  simulation/                        simulate_three_set.py, plot_simulation.py, plot_simulation_supplement.R
  audit_inputs.py, make_analysis_cohort.py, diagnose_exposure.py
  preprocessing/                     TCGA input construction, RNA normalization, missingness report
results/  tmp/                       created at run time (not tracked)
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

- Python ≥ 3.9 with numpy, pandas, scipy, torch, torch-geometric, matplotlib (and scikit-learn for diagnostics)
- R ≥ 4.1 with dplyr, tidyr, readr, ggplot2, patchwork, jsonlite, ggpubr
- CPU is sufficient; memory use grows with patients × edges (TCGA run: about 4 GB of output for 50 fits)

## Citation

If you use ExposureGAT, please cite: *[citation to be added]*

## License

*[license to be added]*
