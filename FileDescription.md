# ExposureGAT file description

This document describes the input files ExposureGAT expects and every output file it writes.
All outputs of one run are written to `results/<name>/`; intermediate and downloaded files go to `tmp/`.

**Naming convention in output columns.** Column names use *reference* for the normal (`NT`) sample and
*transition* for the tumor (`TP`) sample of each pair; *transition* refers to the tumor-minus-normal
difference measured at a single time point, not longitudinal progression. *Exposure* refers to the
column given by `--exposure_col` (default `smoking_intensity`).

---

## 1. Inputs

One folder per dataset. The demo input (`bash exposureGAT.sh demo`) is written to `tmp/demo_data/` and
can be used as a template; the TCGA input is downloaded to `tmp/tcga_data/`.

```text
input_folder/
├── X_rna_gene_features.csv      gene expression        (samples × genes)
├── X_meth_gene_features.csv     promoter methylation   (samples × genes)
├── X_cnv_gene_features.csv      copy number            (samples × genes)
├── X_cpg_features.csv           CpG methylation        (samples × CpGs)
├── edges_ppi.csv                gene–gene edges
├── edges_cpg_to_gene.csv        CpG → gene edges
├── map_gene_nodes.csv           gene identifier → node index
├── map_cpg_nodes.csv            CpG identifier → node index
├── metadata_samples.csv         sample, patient, tissue and exposure information
└── sim_ground_truth_nodes.csv   simulated data only
```

### Node indexing

All nodes share one index space. Gene nodes come first (indices `0 … G−1`, in the column order of the
gene feature matrices) and CpG nodes follow (indices `G … G+C−1`, in the column order of
`X_cpg_features.csv`). Edge files refer to nodes by these indices.

### Feature matrices

All four matrices share the same layout:

- **Rows:** samples. The first column (unnamed) holds the sample identifier, matching `sample_id` in `metadata_samples.csv`.
- **Row order:** identical in all four files.
- **Columns:** one per node.
- **Missing values:** not allowed in RNA and methylation; missing CNV values are read as 0 (diploid) and missing CpG values as 0.5.

Each patient contributes two rows: one normal (reference) sample and one tumor sample. Features are scaled
inside the model (robust median/IQR per gene, fitted on training patients only, followed by tanh), so
supply values unscaled, in the units described below.

#### `X_rna_gene_features.csv`

| Column | Meaning |
|---|---|
| first column | Sample identifier |
| one column per gene (e.g. `1543`) | Gene expression, log scale. TCGA: TMM-normalized log2 counts per million (edgeR, prior count 1). |

Gene column names are the gene identifiers (TCGA: Entrez Gene IDs). They define the gene nodes and their
order, and must be identical, in the same order, in the methylation and CNV matrices.

#### `X_meth_gene_features.csv`

| Column | Meaning |
|---|---|
| first column | Sample identifier |
| one column per gene | Mean beta value (0–1) of the probes in the gene promoter (450K annotation groups TSS1500, TSS200, 5′UTR, first exon). Genes without promoter probes: within-sample median. |

#### `X_cnv_gene_features.csv`

| Column | Meaning |
|---|---|
| first column | Sample identifier |
| one column per gene | Copy number centered at diploid and compressed: sign(x)·log2(\|x\|+1), where x = absolute copy number − 2. 0 = diploid, +1 = one-copy gain, −1 = one-copy loss, −1.58 = homozygous deletion. |

#### `X_cpg_features.csv`

| Column | Meaning |
|---|---|
| first column | Sample identifier |
| one column per CpG (e.g. `cg05575921`) | CpG methylation beta value (0–1). Used without scaling. |

### Edge files

Edges are directed. Edge weights are read for completeness but **not used** by the model: all edges enter
with weight 1, and patient-specific edge importance is learned through attention.

#### `edges_ppi.csv`

| Column | Meaning |
|---|---|
| `from_idx` | Source gene-node index (0 … G−1) |
| `to_idx` | Target gene-node index (0 … G−1) |
| `weight` | Edge weight (not used; TCGA: 1.0) |

Undirected interactions (TCGA: STRING v12, combined score ≥ 700) must be listed in both directions, one row each.

#### `edges_cpg_to_gene.csv`

| Column | Meaning |
|---|---|
| `from_idx` | CpG node index. Either the global index (G … G+C−1) or the CpG's position in `X_cpg_features.csv` (0 … C−1); the loader adds G automatically in the second case. |
| `to_idx` | Target gene-node index (0 … G−1) |
| `weight` | Edge weight (not used; TCGA: 1.0) |
| `cpg_id` | CpG identifier (reference only) |
| `entrez_id` | Identifier of the target gene (reference only) |

List each CpG → gene link once; the loader adds the gene → CpG direction. TCGA: each CpG is linked to the
gene with the nearest transcription start site (hg38).

### Node maps

These files document the node indices. The model takes node identity from the feature-matrix columns, so the
maps are not required for fitting, but they must agree with the edge files if supplied.

#### `map_gene_nodes.csv`

| Column | Meaning |
|---|---|
| `entrez_id` | Gene identifier (same as the gene column name in the feature matrices) |
| `node_idx` | Node index (0 … G−1) |
| `node_type` | 1 = gene |

#### `map_cpg_nodes.csv`

| Column | Meaning |
|---|---|
| `cpg_id` | CpG identifier (same as the column name in `X_cpg_features.csv`) |
| `node_idx` | Node index (G … G+C−1) |
| `node_type` | 0 = CpG |

### `metadata_samples.csv`

One row per sample (two rows per patient).

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Sample identifier matching the feature-matrix rows (`sample_barcode` is also accepted) |
| `patient_id` | yes | Patient identifier shared by the two samples of a pair |
| `tissue_type` | yes | `NT` for the normal (reference) sample, `TP` for the tumor sample. Any value containing NT or TP is accepted. |
| `smoking_intensity` | yes (default exposure) | Exposure value, identical in both rows of a patient. TCGA: pack-years. Another column can be used with `--exposure_col`. |
| `smoking_status` | recommended | Recorded smoking status (e.g. `Lifelong Non-Smoker`, `Current Smoker`). Used to exclude zero/missing exposure without a never-smoker status, and in cohort tables. |
| `age`, `gender` | optional | Cohort characteristics tables |
| `sample_index`, `smoking_rank`, `vital_status`, `diagnosis` | optional | Written by the TCGA build script: row index, ordinal smoking status (1 never … 4 current), vital status, 1 = tumor / 0 = normal. Not used by the model. |

Rules applied when the data are loaded:

- Samples not present in the feature matrices are ignored.
- Rows with missing exposure are dropped. Missing exposure must **not** be coded as zero; use
  `src/make_analysis_cohort.py` to exclude patients whose zero exposure is not confirmed (e.g. not a
  lifelong non-smoker).
- Only patients with both an `NT` and a `TP` sample are analyzed. If a patient has several samples of one
  type, the first by sample identifier is used.
- Patient-level columns (exposure, smoking status, age, gender) must be identical in both rows of a patient.
- At least 10 patients with varying exposure are required.

### `sim_ground_truth_nodes.csv` (simulated data only)

| Column | Meaning |
|---|---|
| `node_idx` | Node index |
| `node_label` | `GENE:<id>` or `CpG:<id>` |
| `set_membership` | `A` = exposure-linked (tumor–normal change scales with exposure), `C` = fixed tumor shift unrelated to exposure, `B` = background |

When present, recovery metrics (AUROC of A vs C and related figures) are produced.

---

## 2. Run settings and logs: `logs/`

| File | Content |
|---|---|
| `settings.txt` | Command, input folder and all settings used. |
| `1_audit.log` … `5_reproducibility.log` | Console output of each pipeline step. |

---

## 3. Model outputs: `model/`

Written by `exposureGAT.py`. A run consists of initialization fits (all patients, different seeds),
subset fits (random 80% of patients) and permutation fits (shuffled exposure).
*Priority* = absolute contribution to the patient's fitted exposure score. *Rank* 1 = highest priority.

### 3.1 Run description

| File | Content |
|---|---|
| `config.json` | All model settings. |
| `runs.csv` | One row per fit: `kind` (initialization / subset / permutation), `repeat`, `seed`, `n_fit` (patients used). |
| `reference_model.pt` | Fitted parameters, cohort reference and scalers of the first initialization fit. |
| `extended_score_metadata.json` | Definitions of edge scores and the edge top-fraction used. |
| `edge_index.csv` | One row per directed edge: `edge_id, source_idx, target_idx, source_label, target_label, edge_type` (`gene_to_gene`, `cpg_to_gene`, `gene_to_cpg`). |
| `simulation_recovery.json` | Simulated data only: AUROC of cohort node priority for set A vs set C. |

### 3.2 Per-fit folders: `initialization_k/`, `subset_k/`, `permutation_k/`

| File | Content |
|---|---|
| `history.csv` | Per epoch: total `loss`, `exposure_loss`, `transition_loss`. |
| `fit_membership.csv` | Per patient: `included_in_fit`, `supplied_exposure` (shuffled in permutation fits). |
| `patient_scores.csv` | Per patient: fitted exposure-rank score, transition score, total absolute contribution, input change magnitude. |
| `attributions.npz` | Patient × node arrays: `score` (priority), `signed` (signed contribution), `neighbor_difference`, `patient_score`. |
| `edge_attributions.npy` | Patient × edge structured array of all edge fields (see 3.4). |
| `node_transition.npy` | Patient × node paired embedding-change magnitude. |

### 3.3 Node scores

| File | Content |
|---|---|
| `patient_node_scores.csv` | One row per patient × node. Columns below. |
| `patient_top_nodes.csv` | Same columns, the `top_k` highest-priority nodes per patient. |
| `cohort_summary.csv` | One row per node: `median_patient_rank` (cohort ranking), `mean_initialization_stability`, `mean_subset_stability`, `mean_null_fraction`, `topk_patient_prevalence` (fraction of patients with the node in their top list). |
| `patient_node_transition_scores.csv.gz` | Per patient × node: paired embedding-change magnitude (`transition_score`). |
| `cohort_node_transition_scores.csv` | Per node: median and mean transition score across patients. |

Columns of `patient_node_scores.csv`:

| Column | Meaning |
|---|---|
| `median_priority` | Median priority across initialization fits. |
| `median_rank` | Median within-patient rank across initialization fits. |
| `initialization_topk_frequency` | Fraction of initialization fits with the node in the patient's top-k. |
| `subset_topk_frequency` | Same across subset fits. |
| `subset_included_frequency` / `subset_excluded_frequency` | Same, restricted to subset fits that included / held out the patient. |
| `null_fraction_at_least_as_high_rank` | Fraction of permutation fits ranking the node at least as high as the consensus rank. Descriptive; not a p-value. |

### 3.4 Edge scores

| File | Content |
|---|---|
| `patient_edge_scores.csv.gz` | One row per patient × directed edge. Columns below. |
| `patient_top_edges.csv.gz` | Same columns, the top fraction of edges per patient (`edge_top_fraction`, default 10%). |
| `cohort_edge_scores.csv` | One row per edge: `median_patient_rank`, `cohort_exposure_priority` (= −median rank), `median_patient_priority`, `median_patient_message_transition_norm`, stability and null summaries, `topk_patient_prevalence`. |
| `patient_node_edge_summary.csv.gz` | Per patient × node: number of selected incoming/outgoing edges and summed incoming/outgoing edge priority. |
| `cohort_node_edge_summary.csv` | The same averaged over patients. |

Columns of `patient_edge_scores.csv.gz` (medians across initialization fits):

| Column | Meaning |
|---|---|
| `median_priority` | Edge priority: absolute change in the patient's fitted score when the edge's message is suppressed (all parameters fixed). |
| `median_signed_exposure_contribution` | The same change with sign. |
| `median_message_transition_norm` | Edge delta: size of the change in the message along the edge between normal and tumor. |
| `median_reference_attention`, `median_transition_attention` | Attention weight in normal and tumor. |
| `median_attention_change`, `median_gated_attention_change` | Tumor minus normal attention (without / with edge gate). |
| `median_rank`, stability and null columns | As for nodes. |

Edge priorities and edge deltas describe how existing connections contribute in each patient; they do not
indicate gained or lost interactions.

### 3.5 Patient scores

| File | Content |
|---|---|
| `patient_scores.csv` | Per patient, summarized over initialization fits (median, min, max, SD): `predicted_exposure_rank_score` (in-sample fit, standardized rank units, not pack-years), `transition_score`, `exposure_contribution_magnitude`, `scaled_input_transition_rms`. |
| `patient_scores_by_run.csv` | The same per fit, with `kind`, `repeat`, `included_in_fit`, `supplied_exposure`. |

---

## 4. Summary tables: `csv/`

### 4.1 `csv/input_audit/`

| File | Content |
|---|---|
| `audit_summary.json` | Patient and sample counts, number with zero exposure by smoking status, unavailable covariates. |
| `exported_patient_audit.csv` | Per patient metadata with a zero-exposure flag. |
| `smoking_status_by_zero_exposure.csv` | Cross-table of smoking status and zero exposure. |
| `exported_molecular_audit.csv` | Per modality: samples, features, non-finite values, constant rows. |
| `cnv_profile_flags.csv`, `cnv_flags_by_tissue.csv` | Samples with all-zero (all-diploid) or constant copy-number profiles. |

### 4.2 `csv/exposure_associated/`

Node priority measures contribution to the fitted score, not association with exposure. These tables
identify nodes whose contributions track exposure beyond what shuffled-label fits produce.

| File | Content |
|---|---|
| `Exposure_associated_node_scores_<E>epochs.csv` | Per node: `normalized_priority`, `exposure_rho` (Spearman of signed contribution with exposure, median over initializations; `_min`/`_max` across fits), `null_percentile` (share of pooled shuffled-label alignments exceeded), `exposure_associated_score` (= normalized priority × null percentile), `exposure_associated` (null percentile ≥ threshold, default 0.99), `input_delta_rho` (per-feature Spearman of the paired input difference with exposure, strongest modality). |
| `Exposure_associated_patient_scores_<E>epochs.csv` | Per patient × exposure-associated node: priority, signed contribution, patient-level exposure-associated score. |
| `Exposure_associated_summary.csv` | Per node type: number passing, number expected by chance, median number passing per real and per shuffled fit. |
| `Exposure_associated_overlap_checks.csv` | Overlap of the exposure-associated set with the top-priority list and with the top per-feature correlation list, and the overlap expected by chance. |
| `Exposure_associated_notes.txt` | Definitions and threshold. |

### 4.3 `csv/reproducibility/`

| File | Content |
|---|---|
| `Held_out_summary.csv` | Fitted and held-out agreement with exposure, top-node retention across seeds and splits (patient included vs held out), chance retention, cohort ranking agreement. |
| `Fit_by_seed.csv` | In-sample Spearman per initialization fit. |
| `Held_out_by_split.csv` | Per subset fit: training and held-out patient counts and Spearman. |
| `Held_out_pooled_by_patient.csv` | Per patient: mean score over fits that held the patient out. |
| `Patient_top_node_stability.csv` | Per patient: top-node retention across seeds, splits with the patient included, splits with the patient held out. |
| `Cohort_ranking_stability.csv` | Rank agreement and top-k overlap for every pair of seeds and every pair of splits. |
| `Held_out_response_text.txt` | Text summary with values filled in. |

### 4.4 `csv/figures/`

Source data for the figures and cohort tables (file names starting with `TCGA_` are produced for any dataset).

| File | Content |
|---|---|
| `patient_plot_data.csv`, `node_plot_data.csv`, `edge_plot_data.csv` | Data behind the patient, node and edge figures. |
| `Fit_check_summary.csv`, `final_fit_correlations.csv` | Fitted agreement with real vs shuffled labels per fit. |
| `Node_agreement_summary.csv`, `Node_agreement_real_vs_shuffled_*.csv` | Node ranking agreement between real-label fits vs between real- and shuffled-label fits (patient and cohort level). |
| `patient_correlations.csv` | Per dataset: Spearman of transition score and fitted score with exposure. |
| `TCGA_available_cohort_counts.csv`, `TCGA_available_cohort_numeric.csv`, `TCGA_supplied_patient_metadata.csv` | Cohort characteristics from the supplied metadata. |
| `TCGA_exposure_group_definitions.csv`, `TCGA_patient_plot_data.csv` | Exposure tertile groups (ties kept together) and per-patient data. |
| `TCGA_top10percent_genes.csv`, `Top_cohort_gene_nodes_all_ranked.csv`, `Top_cohort_CpG_nodes_all_ranked.csv` | Cohort gene and CpG rankings (top 10% and full list). |
| `TCGA_top_edges.csv`, `TCGA_top_CpG_gene_edges.csv`, `TCGA_top_CpG_nodes.csv` | Top cohort edges, CpG–gene edges and CpG nodes. |
| `TCGA_selection_boundaries.csv`, `TCGA_Supp7_cutoffs.csv` | Top-fraction cutoffs used for selections. |
| `TCGA_patient_node_edge_summaries.csv`, `TCGA_Supp8_patient_values.csv`, `TCGA_node_edge_correlation.csv` | Per-patient node and edge summaries and their relationship. |
| `TCGA_Supp10_*` | Values behind the CpG heatmap. |
| `edge_selection_vs_graph_baseline.csv`, `actual_input_degrees.csv`, `node_edge_summary_plot_data.csv`, `edge_attention_plot_data.csv`, `patient_edge_heatmap_data.csv` | Edge selection by endpoint class, input degrees, connectivity and attention summaries. |
| `observed_magnitude_data.csv`, `magnitude_exposure_correlations.csv` | Raw paired input change per patient and its correlation with exposure. |
| `AUROC_by_dataset.csv`, `AUROC_summary.csv`, `top_node_recovery.csv` | Simulated data only: recovery of planted node sets. |
| `analysis_notes.txt`, `TCGA_limitations.txt`, `R_session_info.txt` | Scope notes and R session information. |

---

## 5. Figures: `plots/`

Each figure is saved as PDF and PNG.

| Folder / file | Content |
|---|---|
| `figures/Fit_check_real_vs_shuffled` | Fitted vs recorded exposure rank, next to fitted agreement for real and shuffled labels. |
| `figures/Node_agreement_real_vs_shuffled` | Node ranking agreement between real-label fits vs real- and shuffled-label fits (patient and cohort). |
| `figures/Predicted_vs_actual_smoking_rank`, `Predicted_rank_vs_recorded_smoking` | In-sample fitted score vs exposure. |
| `figures/Patient_transition_vs_normalized_exposure` | Patient paired-change magnitude vs exposure. |
| `figures/Figure3B_*`, `Figure3C_*`, `Figure3_ABC_no_pathway_panel` | Patient node contributions and edge strength by exposure group (descriptive). |
| `figures/Top_cohort_gene_nodes`, `Top_cohort_CpG_nodes`, `Additional_patient_node_heatmap`, `Additional_patient_edge_heatmap` | Top cohort nodes and patient-level node/edge priority heatmaps. |
| `figures/Additional_node_confidence`, `Additional_edge_confidence`, `Additional_final_fit_correlations` | Stability and null summaries; per-fit agreement. |
| `figures/Supp1_training_losses`, `TCGA_Supp6_training` | Training loss curves. |
| `figures/Supp3_*`, `Supp4_*`, `Supp5_*`, `Edge_message_change_vs_priority`, `Additional_attention_change_vs_edge_priority` | Node transition scores, edge selection and priority, connectivity, edge delta and attention vs priority. |
| `figures/TCGA_Supp7_*`, `TCGA_Supp8_*`, `TCGA_Supp10_*` | Score distributions and cutoffs, node–edge relationship, CpG summaries. |
| `figures/Reviewer_magnitude_and_connectivity` | Raw paired change vs exposure; input connectivity. |
| `figures/Figure2_combined`, `Reviewer_A_vs_C_recovery`, `Reviewer_node_recall`, `Reviewer_top_nodes_vs_expectation`, `Supp2_*` | Simulated data only: recovery of planted node sets. |
| `exposure_associated/Exposure_associated_genes_heatmap`, `..._CpGs_heatmap` | Signed contributions of exposure-associated nodes per patient, ordered by exposure. |
| `exposure_associated/Exposure_associated_counts_real_vs_shuffled` | Number of nodes passing the exposure-alignment rule in each real and shuffled fit vs chance. |
| `reproducibility/Held_out_and_stability` | Fitted and held-out agreement; top-node retention across seeds and splits; cohort ranking agreement. |

---

## 6. Other scripts and their outputs

| Script | Output |
|---|---|
| `src/simulation/simulate_data.py` | A simulated input folder (Section 1), used by `exposureGAT.sh demo`. |
| `src/make_analysis_cohort.py` | A copy of an input folder with excluded patients removed, plus `cohort_selection.csv`. |
| `src/diagnose_exposure.py` | Training-length (`--stage opt`), loss-weight (`loss`), capacity (`capacity`) and input-signal (`ceiling`) diagnostics: `trajectory.csv`, `checkpoint_summary.csv`, `ceiling_cv_ridge.csv`, figures. |
| `src/analysis/node_heatmaps.R` | Heatmaps of top-priority genes and CpGs across several runs (one page per run). |
| `src/preprocessing/build_tcga_inputs.R` | TCGA inputs (Section 1) from GDC data. |
| `src/preprocessing/normalize_rna.R` | Input folder with TMM log2-CPM RNA, plus `rna_normalization_checks/`. |
| `src/preprocessing/missingness_report.R` | `Missingness_table.csv`: replaced methylation, copy-number and CpG values, and exposure zeros by smoking status. |
