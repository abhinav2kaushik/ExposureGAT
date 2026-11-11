# ExposureGAT file description

This document describes the input files ExposureGAT expects and every output file it writes.
All outputs of one run are written to `results/<name>/`; intermediate and downloaded files go to `tmp/`.

**Naming convention in output columns.** Column names use *reference* for the normal (`NT`) sample and
*transition* for the tumor (`TP`) sample of each pair; *transition* refers to the tumor-minus-normal
difference measured at a single time point, not longitudinal progression. *Exposure* refers to the
column given by `--exposure_col` (default `smoking_intensity`).

---

## 1. Inputs

One folder per dataset. The demo input (`python3 exposureGAT.py demo`) is written to `tmp/demo_data/` and
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

Written by `src/exposuregat_model.py` (called by `exposureGAT.py`). A run consists of initialization fits (all patients, different seeds),
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

### 4.2 `csv/weighted_node_score/` (main output)

The weighted node score ranks nodes by how consistently their patient contributions align with the supplied
exposure, with a small weight on paired-transition magnitude: W = max(A, 0) × ((1 − γ) + γ·T), where A is the
Spearman correlation of patient contributions with exposure (median over initializations), T the percentile of
the median paired-transition magnitude, and γ = `--gamma` (default 0.10). It is descriptive and within cohort.
It is computed for every run (`demo`, `simulation`, `tcga`, `run`) and for every node type in the graph: each
type is scored and ranked separately (T is the percentile among nodes of the same type), so gene ranks do not
depend on how many CpGs are in the graph. Gene-node files carry no suffix; files for other node types carry the
type as suffix (e.g. `weighted_node_scores_CpG.csv`, with `feature_id` in place of `entrez_id`).

| File | Content |
|---|---|
| `weighted_node_scores.csv` | **Main output.** Per node: `rank`, `node_idx`, `node_label`, `entrez_id`, `weighted_node_score`, `exposure_alignment` (A), `transition_magnitude`, `transition_percentile` (T), `top_fraction` (in the top `--top_fraction`, default 10%). Gene nodes. |
| `weighted_node_scores_CpG.csv` | The same for CpG nodes, ranked among CpGs (when CpG nodes are present). |
| `weighted_node_scores_by_fit.csv.gz`, `weighted_node_scores_by_fit_CpG.csv.gz` | The same quantities for every fit (initializations, patient subsets, permutations), with `kind`, `run` and `selected` (in the top fraction of that fit). |
| `patient_weighted_contributions.csv.gz`, `patient_weighted_contributions_CpG.csv.gz` | Per patient × top node: weighted score × absolute median contribution (patient-level summary). |
| `top_genes_signed_contributions.csv`, `top_CpGs_signed_contributions.csv` | Signed contributions of the top 30 nodes per patient, patients ordered by exposure (Figure 3B data for genes). |
| `top_gene_alignment_bins.csv`, `top_CpG_alignment_bins.csv` | Distribution of exposure alignment among the top 5%, 10% and 20% nodes (Figure 3C data for genes). |
| `settings.json` | γ, top fraction, node types with their number of nodes and top nodes, number of fits. |

### 4.3 `csv/label_controls/`

Real-label fits compared with shuffled-label fits (two-sided run-group randomization, 9,999 draws; Holm
adjustment across the five primary comparisons).

| File | Content |
|---|---|
| `patient_ranking_agreement.csv`, `cohort_ranking_agreement.csv` | Rank Spearman and top-100 Jaccard for real–real and real–shuffled fit pairs. |
| `fixed_top5_correlations.csv` | Agreement of each fit with the fixed top 5% genes of the real-label consensus. |
| `alignment_counts.csv`, `alignment_count_summary.csv` | Number of genes with exposure alignment ρ ≥ 0.30 per real and shuffled fit. |
| `statistical_comparisons.csv` | Effect, p, test and Holm-adjusted p per comparison. |
| `summary.json` | Settings and counts. |

### 4.4 `csv/genes/`

| File | Content |
|---|---|
| `direction_of_change_top_genes.csv` | Top genes by weighted node score: RNA median transition − reference, fraction higher in transition, direction; promoter-CpG methylation (complete pairs, median change, direction). |
| `marker_ranks.csv` | Ranks and scores of the marker genes and CpGs given by `--markers`/`--marker_cpgs` (TCGA: AHRR, CYP1A1, CYP1B1, FXYD1, MYO1G and seven smoking-associated CpGs); genes are ranked among genes and CpGs among CpGs; markers not present in the graph are listed as absent. |

### 4.5 `csv/gene_signatures/` (only when MSigDB GMT files are supplied)

| File | Content |
|---|---|
| `gene_signature_table.csv` | Per signature in `resources/signature_catalog.csv`: collection, size, genes in the graph, genes in the top 10%, expected count, fold representation, retained/excluded with reason. |
| `signature_summary.csv`, `gene_set_union_summary.csv` | Counts for the smoking- and lung cancer-related signature unions. |
| `gene_set_membership.csv` | Signature membership of every gene, with its weighted node score, rank, exposure alignment and top-fraction flag (empty when the gene is not in the graph). |
| `signature_counts_by_fit.csv`, `real_minus_shuffled.csv` | Signature counts in the top genes of each real and shuffled fit. |
| `weighted_score_density_data.csv` | Density curves of weighted node scores of signature genes in the top 5%, 10% and 20% (Supplementary Figure 9 data). |

### 4.6 `csv/reproducibility/`

| File | Content |
|---|---|
| `Held_out_summary.csv` | Fitted and held-out agreement with exposure, top-node retention across seeds and splits (patient included vs held out), chance retention, cohort ranking agreement. |
| `Fit_by_seed.csv` | In-sample Spearman per initialization fit. |
| `Held_out_by_split.csv` | Per subset fit: training and held-out patient counts and Spearman. |
| `Held_out_pooled_by_patient.csv` | Per patient: mean score over fits that held the patient out. |
| `Patient_top_node_stability.csv` | Per patient: top-node retention across seeds, splits with the patient included, splits with the patient held out. |
| `Cohort_ranking_stability.csv` | Rank agreement and top-k overlap for every pair of seeds and every pair of splits. |

### 4.7 `csv/figures/`

Source data for the figures and cohort tables (file names starting with `TCGA_` are produced for any dataset).

| File | Content |
|---|---|
| `patient_plot_data.csv`, `node_plot_data.csv`, `edge_plot_data.csv` | Data behind the patient, node and edge figures. |
| `Fit_check_summary.csv`, `final_fit_correlations.csv` | Fitted agreement with real vs shuffled labels per fit. |
| `Node_agreement_summary.csv`, `Node_agreement_real_vs_shuffled_*.csv` | Node ranking agreement between real-label fits vs between real- and shuffled-label fits (patient and cohort level). |
| `patient_correlations.csv` | Per dataset: Spearman of transition score and fitted score with exposure. |
| `TCGA_available_cohort_counts.csv`, `TCGA_available_cohort_numeric.csv`, `TCGA_supplied_patient_metadata.csv` | Cohort characteristics from the supplied metadata. |
| `TCGA_exposure_group_definitions.csv`, `TCGA_patient_plot_data.csv` | Exposure tertile groups (ties kept together) and per-patient data. |
| `TCGA_top10percent_genes.csv`, `Top_cohort_gene_nodes_all_ranked.csv`, `Top_cohort_CpG_nodes_all_ranked.csv` | Cohort gene and CpG rankings by node priority (top 10% and full list); the main ranking is in `csv/weighted_node_score/`. |
| `TCGA_top_edges.csv`, `TCGA_top_CpG_gene_edges.csv`, `TCGA_top_CpG_nodes.csv` | Top cohort edges, CpG–gene edges and CpG nodes. |
| `TCGA_selection_boundaries.csv`, `TCGA_Supp7_cutoffs.csv` | Top-fraction cutoffs used for selections. |
| `TCGA_patient_node_edge_summaries.csv`, `TCGA_Supp8_patient_values.csv`, `TCGA_node_edge_correlation.csv` | Per-patient node and edge summaries and their relationship. |
| `TCGA_Supp10_*` | Values behind the CpG heatmap (produced when CpG nodes are present). |
| `edge_selection_vs_graph_baseline.csv`, `actual_input_degrees.csv`, `node_edge_summary_plot_data.csv`, `edge_attention_plot_data.csv`, `patient_edge_heatmap_data.csv` | Edge selection by endpoint class, input degrees, connectivity and attention summaries. |
| `observed_magnitude_data.csv`, `magnitude_exposure_correlations.csv` | Raw paired input change per patient and its correlation with exposure. |
| `AUROC_by_dataset.csv`, `AUROC_summary.csv`, `top_node_recovery.csv` | Simulated data only: recovery of planted node sets. |
| `analysis_notes.txt`, `TCGA_limitations.txt`, `R_session_info.txt` | Scope notes and R session information. |

---

## 5. Figures: `plots/`

Each figure is saved as PDF and PNG.

| Folder / file | Content |
|---|---|
| `weighted_node_score/Top_genes_signed_contributions_heatmap` | Signed contributions of the top 30 genes by weighted node score per patient, ordered by exposure (Figure 3B). |
| `weighted_node_score/Top_CpGs_signed_contributions_heatmap`, `Top_CpG_exposure_alignment` | The same for CpG nodes (when present). |
| `weighted_node_score/Top_gene_exposure_alignment` | Exposure alignment of the top 5%, 10% and 20% genes (Figure 3C). |
| `label_controls/Ranking_agreement_real_vs_shuffled` | Ranking agreement between real-label fits vs real- and shuffled-label fits (Supplementary Figure 10). |
| `label_controls/Exposure_aligned_gene_counts`, `Fixed_top5_agreement_with_shuffled_fits` | Aligned-gene counts and fixed top 5% agreement, real vs shuffled (Supplementary Figure 11). |
| `gene_signatures/Smoking_signatures`, `Lung_cancer_signatures` | Representation of signature genes in the top 10% (Supplementary Figures 12–13). |
| `gene_signatures/Weighted_score_density_*` | Weighted scores of signature genes vs all genes (Supplementary Figure 9). |
| `figures/Fitted_vs_recorded_exposure_rank`, `Fitted_rank_vs_recorded_exposure` | In-sample fitted score vs exposure (Figure 3A). |
| `figures/Fit_check_real_vs_shuffled`, `Node_agreement_real_vs_shuffled` | Fitted agreement and node ranking agreement, real vs shuffled labels. |
| `figures/Patient_transition_vs_normalized_exposure` | Patient paired-transition magnitude vs exposure. |
| `figures/Patient_node_contributions_by_exposure_group`, `Patient_edge_strength_by_exposure_group`, `Patient_summaries_combined` | Patient node and edge summaries by exposure tertile (descriptive). |
| `figures/Top_cohort_gene_nodes`, `Additional_patient_node_heatmap`, `Additional_patient_edge_heatmap` | Top cohort genes by node priority; patient node/edge priority heatmaps. |
| `figures/Additional_node_confidence`, `Additional_edge_confidence`, `Additional_final_fit_correlations` | Stability summaries; per-fit agreement. |
| `figures/Supp1_training_losses`, `TCGA_Supp6_training` | Training loss curves (Supplementary Figure 6 for TCGA). |
| `figures/Supp2_*`, `Supp3_*`, `Supp4_*`, `Supp5_*`, `Edge_message_change_vs_priority`, `Additional_attention_change_vs_edge_priority` | Node priority, paired-transition scores, edge selection and priority, connectivity, message and attention change vs edge priority. |
| `figures/TCGA_Supp7_*`, `TCGA_Supp8_*` | Node and edge priority distributions and cutoffs; node–edge relationship (Supplementary Figures 7–8). |
| `figures/Figure2_combined`, `Simulation_*` | Simulated data only: recovery of planted node sets. |
| `reproducibility/Held_out_and_stability` | Fitted and held-out agreement; top-node retention across seeds and splits; cohort ranking agreement. |

### 5.1 Simulation benchmark: `results/simulation/` (`python3 exposureGAT.py simulation`)

| File | Content |
|---|---|
| `manifest.csv` | Simulated datasets (one per seed) and their result folders `seed_<n>/`. |
| `seed_<n>/model/` | Model outputs per dataset (Section 3). |
| `seed_<n>/weighted_node_score/` | Weighted node scores per dataset (Section 4.2). |
| `figures/Figure_2` | AUROC of exposure-linked (A) vs fixed-shift (C) nodes, top-k recovery and patient agreement (Figure 2; curves averaged across datasets). |
| `figures/AUROC_by_dataset.csv`, `AUROC_summary.csv`, `repetition_AUROC.csv`, `topk_recovery.csv`, `node_scores.csv`, `top100_seed_<n>.csv` | Data behind Figure 2. |
| `figures/Supplement_1_training` … `Supplement_5_robustness`, `design_audit.csv`, `patient_magnitude_audit.csv`, `actual_node_degrees.csv` | Training curves, node priority by class, cohort node priority by class, design checks and robustness (Supplementary Figures 1–5). |
| `supplementary_figures/` | Additional simulation figures from `plot_simulation_supplement.R`. |

---

## 6. Other scripts and their outputs

| Script | Output |
|---|---|
| `src/simulation/simulate_three_set.py` | A simulated input folder (Section 1) with planted exposure-linked (A), fixed-shift (C) and background (B) nodes; used by `demo` and `simulation`. |
| `src/simulation/plot_simulation.py`, `plot_simulation_supplement.R` | Figure 2 and simulation supplementary figures (Section 5.1). |
| `src/scoring/weighted_node_score.py` | Section 4.2 from a `model/` folder. |
| `src/analysis/label_controls.py`, `src/plotting/plot_label_controls.R` | Section 4.3 and its figures. |
| `src/analysis/direction_of_change.py`, `marker_ranks.py` | Section 4.4. |
| `src/analysis/gene_signatures.py`, `src/plotting/plot_gene_signatures.R` | Section 4.5 and its figures (needs MSigDB Entrez GMT files). |
| `src/make_analysis_cohort.py` | A copy of an input folder with excluded patients removed, plus `cohort_selection.csv`. |
| `src/diagnose_exposure.py` | Training-length (`--stage opt`), loss-weight (`loss`), capacity (`capacity`) and input-signal (`ceiling`) diagnostics: `trajectory.csv`, `checkpoint_summary.csv`, `ceiling_cv_ridge.csv`, figures. |
| `src/preprocessing/build_tcga_inputs.R` | TCGA inputs (Section 1) from GDC data. |
| `src/preprocessing/normalize_rna.R` | Input folder with TMM log2-CPM RNA, plus `rna_normalization_checks/`. |
| `src/preprocessing/missingness_report.R` | `Missingness_table.csv`: replaced methylation, copy-number and CpG values, and exposure zeros by smoking status. |
