# =============================================================================
# ExposureGAT TCGA: ALL plots in one script
# =============================================================================
# Combines the former plot_all_results.R (shared node/edge/patient plots) and
# plot_tcga.R (TCGA Figure 3, Supplements 6-10), plus new sections:
#   N1. Fit check figure: in-sample fit next to shuffled-label fits
#   N2. Real-vs-shuffled node agreement (initialization vs permutation runs)
#   N3. Epoch robustness: compare node/edge rankings with a second analysis run
#   N4. Training-length / ceiling / marker diagnostics from diagnose_exposure.py
# No existing plot was removed. No fitting, p-values or model-score changes.
#
# Run from Terminal (via plot_tcga.sh or directly):
#   Rscript plot_tcga.R INPUT_DATA ANALYSIS_DIR PLOT_DIR [DIAGNOSTICS_DIR] [COMPARE_ANALYSIS_DIR]
# Use '-' to skip an optional folder. Example:
#   Rscript plot_tcga.R "$C" "$V/TCGA_clean59_e600/analysis" "$V/TCGA_clean59_e600/plots" \
#     "$V/diagnostics" "$V/TCGA_clean59_e400/analysis"
# INPUT_DATA must be the folder the model was fitted on (e.g. GNN_Input_Data_clean59).
# Optional: MANIFEST=/path/manifest.csv (columns seed,data_dir,analysis_dir) for
# multiple simulated datasets, as before.
# install.packages(c('dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite'))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)

# =============================================================================
# Settings: edit these paths when running interactively in RStudio
# =============================================================================
data_dir <- '/path/to/GNN_Input_Data_clean59'
analysis_dir <- '/path/to/TCGA_clean59_e600/analysis'
output_dir <- '/path/to/TCGA_clean59_e600/plots'
diagnostics_dir <- ''       # folder holding diagnose_exposure.py outputs, or ''
compare_analysis_dir <- ''  # second analysis folder (e.g. 400 epochs), or ''
manifest_file <- Sys.getenv('MANIFEST', '')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  data_dir <- args[1]
  analysis_dir <- args[2]
  output_dir <- args[3]
}
if (length(args) >= 4 && args[4] != '-') diagnostics_dir <- args[4]
if (length(args) >= 5 && args[5] != '-') compare_analysis_dir <- args[5]
if (!dir.exists(data_dir)) stop('INPUT_DATA folder not found: ', data_dir)
if (!dir.exists(analysis_dir)) stop('ANALYSIS_DIR folder not found: ', analysis_dir)
if (nzchar(diagnostics_dir) && !dir.exists(diagnostics_dir)) stop('DIAGNOSTICS_DIR not found: ', diagnostics_dir)
if (nzchar(compare_analysis_dir) && !dir.exists(compare_analysis_dir)) stop('COMPARE_ANALYSIS_DIR not found: ', compare_analysis_dir)

# #############################################################################
# PART 1: shared node, edge and patient plots (formerly plot_all_results.R)
# #############################################################################
cohort_top_k <- 50
heatmap_top_nodes <- 40
heatmap_top_edges <- 20
plot_patient_edges <- TRUE # Set FALSE to skip reading the large patient-edge table.
exposure_column <- 'smoking_intensity'
class_colors <- c(A = '#ED7658', C = '#8361BB', B = '#939D9D')
kind_colors <- c(initialization = '#2886B5', subset = '#C17A26', permutation = '#A14792')
set.seed(42) # Only controls plotting jitter, never model scores.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Small helpers keep the plotting sections readable.
read_table <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)
save_plot <- function(plot, name, width = 8, height = 5) {
  ggsave(file.path(output_dir, paste0(name, '.png')), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(name, '.pdf')), plot, width = width, height = height)
}
rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = 'spearman')
}
auc <- function(group, score, positive, negative) {
  keep <- group %in% c(positive, negative) & is.finite(score)
  y <- group[keep] == positive
  if (!any(y) || !any(!y)) return(NA_real_)
  (sum(rank(score[keep])[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'),
        legend.position = 'bottom')
theme_set(plot_theme)

# =============================================================================
# Input datasets: one real cohort or multiple independent simulations
# =============================================================================
if (nzchar(manifest_file)) {
  datasets <- read_table(manifest_file) %>% mutate(seed = as.character(seed))
} else {
  datasets <- tibble(seed = '1', data_dir = data_dir, analysis_dir = analysis_dir)
}
stopifnot(all(c('seed', 'data_dir', 'analysis_dir') %in% names(datasets)))
if (anyDuplicated(datasets$seed)) stop('Dataset identifiers must be unique.')
# Manifest first row supplies illustrative panels; never select the best seed.
example_seed <- datasets$seed[1]
node_list <- list(); patient_list <- list(); edge_list <- list()
patient_node_list <- list(); degree_list <- list(); history_list <- list()
run_patient_list <- list(); input_degree_list <- list(); audit_list <- list()
latent_list <- list()

for (i in seq_len(nrow(datasets))) {
  d <- datasets$data_dir[i]; a <- datasets$analysis_dir[i]; id <- datasets$seed[i]
  required <- c('cohort_summary.csv', 'patient_scores.csv', 'cohort_edge_scores.csv',
                'cohort_node_transition_scores.csv', 'cohort_node_edge_summary.csv',
                'patient_node_scores.csv', 'patient_scores_by_run.csv', 'runs.csv')
  missing <- required[!file.exists(file.path(a, required))]
  if (length(missing)) stop(paste('Missing extension outputs in', a, ':', paste(missing, collapse = ', ')))
  config <- fromJSON(file.path(a, 'config.json'))
  exposure_column <- config$exposure_col
  meta <- read_table(file.path(d, 'metadata_samples.csv')) %>%
    transmute(patient_id = as.character(patient_id), actual_exposure = .data[[exposure_column]]) %>% distinct()
  if (anyDuplicated(meta$patient_id)) stop('Different exposure values for the same paired patient.')

  # Average ranks preserve ties. Standardization uses population SD, matching Python.
  p <- read_table(file.path(a, 'patient_scores.csv')) %>%
    mutate(patient_id = as.character(patient_id)) %>% left_join(meta, by = 'patient_id')
  if (any(!is.finite(p$actual_exposure))) stop('Missing/nonfinite exposure after joining patient identifiers.')
  p <- p %>% mutate(seed = id,
    normalized_exposure = (rank(actual_exposure, ties.method = 'average') - 0.5) / n(),
    actual_exposure_rank_z = (normalized_exposure - mean(normalized_exposure)) /
      sqrt(mean((normalized_exposure - mean(normalized_exposure))^2)),
    predicted_exposure_rank_z = predicted_exposure_rank_score_median,
    patient_transition_score = transition_score_median)
  if (any(!is.finite(p$actual_exposure_rank_z))) stop('Exposure must vary within a dataset.')
  patient_list[[i]] <- p

  gt_file <- file.path(d, 'sim_ground_truth_nodes.csv')
  if (file.exists(gt_file)) {
    truth <- read_table(gt_file) %>% select(node_idx, set_membership)
  } else {
    truth <- read_table(file.path(a, 'cohort_summary.csv')) %>%
      transmute(node_idx, set_membership = 'Unlabelled')
  }
  nodes <- read_table(file.path(a, 'cohort_summary.csv')) %>%
    left_join(truth, by = 'node_idx') %>%
    left_join(read_table(file.path(a, 'cohort_node_transition_scores.csv')),
              by = c('node_idx', 'node_label')) %>%
    mutate(seed = id, cohort_priority = -median_patient_rank,
           normalized_priority = (n() - median_patient_rank) / max(n() - 1, 1))
  if (anyNA(nodes$set_membership)) stop('Incomplete node truth mapping.')
  node_list[[i]] <- nodes
  patient_node_list[[i]] <- read_table(file.path(a, 'patient_node_scores.csv')) %>%
    left_join(truth, by = 'node_idx') %>% mutate(seed = id)
  degree_list[[i]] <- read_table(file.path(a, 'cohort_node_edge_summary.csv')) %>%
    left_join(truth, by = 'node_idx') %>% mutate(seed = id)

  # Edge classes describe their endpoint NODE classes, not planted edge truth.
  edge_config <- fromJSON(file.path(a, 'extended_score_metadata.json'))
  edges <- read_table(file.path(a, 'cohort_edge_scores.csv')) %>%
    left_join(truth %>% rename(source_idx = node_idx, source_set = set_membership), by = 'source_idx') %>%
    left_join(truth %>% rename(target_idx = node_idx, target_set = set_membership), by = 'target_idx') %>%
    mutate(seed = id, endpoint_group = paste(pmin(source_set, target_set),
           pmax(source_set, target_set), sep = '-'),
           normalized_priority = (n() - median_patient_rank) / max(n() - 1, 1)) %>%
    arrange(median_patient_rank, edge_id) %>%
    mutate(selected = row_number() <= edge_config$edge_top_k)
  edge_list[[i]] <- edges

  # Audit unique undirected input connections, rather than counting both directions.
  undirected <- edges %>% transmute(u = pmin(source_idx, target_idx),
    v = pmax(source_idx, target_idx)) %>% filter(u != v) %>% distinct()
  input_degree <- bind_rows(undirected %>% transmute(node_idx = u),
                           undirected %>% transmute(node_idx = v)) %>% count(node_idx, name = 'degree')
  input_degree_list[[i]] <- truth %>% left_join(input_degree, by = 'node_idx') %>%
    mutate(degree = replace_na(degree, 0L), seed = id)

  runs <- read_table(file.path(a, 'runs.csv'))
  for (j in seq_len(nrow(runs))) {
    folder <- paste0(runs$kind[j], '_', runs[["repeat"]][j])
    h <- read_table(file.path(a, folder, 'history.csv')) %>%
      mutate(seed = id, kind = runs$kind[j], repetition = runs[["repeat"]][j])
    history_list[[length(history_list) + 1]] <- h
  }
  run_patient_list[[i]] <- read_table(file.path(a, 'patient_scores_by_run.csv')) %>% mutate(seed = id)

  # Observed raw modality magnitudes, preserving NT/TP pairing from metadata.
  sample_meta <- read_table(file.path(d, 'metadata_samples.csv')) %>%
    select(sample_id, patient_id, tissue_type) %>% mutate(across(everything(), as.character))
  for (modality in c('rna', 'meth', 'cnv')) {
    x <- read_table(file.path(d, paste0('X_', modality, '_gene_features.csv')))
    names(x)[1] <- 'sample_id'
    normal_ids <- sample_meta %>% filter(tissue_type == 'NT') %>%
      select(patient_id, normal_id = sample_id)
    tumor_ids <- sample_meta %>% filter(tissue_type == 'TP') %>%
      select(patient_id, tumor_id = sample_id)
    pairs <- inner_join(normal_ids, tumor_ids, by = 'patient_id')
    # Raw audit is skipped if this input uses different tissue labels or multiple pairs.
    if (nrow(pairs) == nrow(p) && !anyDuplicated(pairs$patient_id) &&
        all(pairs$normal_id %in% x$sample_id) && all(pairs$tumor_id %in% x$sample_id)) {
      xm <- as.matrix(x[, -1])
      delta <- xm[match(pairs$tumor_id, x$sample_id), , drop = FALSE] -
               xm[match(pairs$normal_id, x$sample_id), , drop = FALSE]
      audit_list[[length(audit_list) + 1]] <- tibble(seed = id,
        patient_id = pairs$patient_id, measure = paste0(modality, '_L2_squared'),
        magnitude = rowSums(delta^2)) %>% left_join(meta, by = 'patient_id')
    } else message('Skipping raw ', modality, ' audit for ', id, ': expected one NT/TP pair per patient.')
  }
  # RMS is monotonically related to total squared change for fixed feature count.
  audit_list[[length(audit_list) + 1]] <- p %>% transmute(seed, patient_id,
    measure = 'loader_transformed_RMS', magnitude = scaled_input_transition_rms_median, actual_exposure)
  checks_file <- file.path(d, 'design_checks.json')
  if (file.exists(checks_file)) {
    checks <- fromJSON(checks_file)
    for (modality in intersect(c('rna', 'meth', 'cnv'), names(checks))) {
      latent_list[[length(latent_list) + 1]] <- tibble(seed = id, modality,
        latent_L1_range = checks[[modality]]$latent_L1_range,
        latent_L2_squared_range = checks[[modality]]$latent_L2_squared_range)
    }
  }
}

nodes <- bind_rows(node_list); patients <- bind_rows(patient_list)
edges <- bind_rows(edge_list); patient_nodes <- bind_rows(patient_node_list)
degrees <- bind_rows(degree_list); histories <- bind_rows(history_list)
run_patients <- bind_rows(run_patient_list); input_degrees <- bind_rows(input_degree_list)
audit <- bind_rows(audit_list)
example_patients <- patients %>% filter(seed == example_seed)
example_nodes <- nodes %>% filter(seed == example_seed)
example_edges <- edges %>% filter(seed == example_seed)
has_truth <- all(nodes$set_membership %in% c('A', 'B', 'C'))
write_csv(patients, file.path(output_dir, 'patient_plot_data.csv'))
write_csv(nodes, file.path(output_dir, 'node_plot_data.csv'))
write_csv(edges, file.path(output_dir, 'edge_plot_data.csv'))
write_csv(input_degrees, file.path(output_dir, 'actual_input_degrees.csv'))

# =============================================================================
# Transition score vs smoking
# =============================================================================
patient_correlations <- patients %>% group_by(seed) %>% summarise(
  n_patients = n(),
  transition_exposure_spearman = rho(patient_transition_score, actual_exposure),
  predicted_actual_spearman = rho(predicted_exposure_rank_z, actual_exposure),
  .groups = 'drop')
write_csv(patient_correlations, file.path(output_dir, 'patient_correlations.csv'))
transition_rho <- rho(example_patients$patient_transition_score, example_patients$normalized_exposure)
p_transition <- ggplot(example_patients, aes(normalized_exposure, patient_transition_score)) +
  geom_point(alpha = 0.7, color = '#337A96') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black') +
  labs(title = 'Patient transition score vs exposure',
       subtitle = paste0('Dataset ', example_seed, '; Spearman rho = ', round(transition_rho, 3)),
       x = 'Normalized exposure rank (ties averaged)', y = 'Mean paired embedding-change magnitude',
       caption = 'Median across full-cohort initializations; line is a descriptive linear trend.')
save_plot(p_transition, 'Patient_transition_vs_normalized_exposure')

# =============================================================================
# Predicted smoking vs actual smoking: matched rank units
# =============================================================================
predicted_rho <- rho(example_patients$predicted_exposure_rank_z, example_patients$actual_exposure_rank_z)
p_predicted <- ggplot(example_patients, aes(actual_exposure_rank_z, predicted_exposure_rank_z)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_point(alpha = 0.7, color = '#8361BB') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black') +
  labs(title = 'Predicted exposure rank vs actual exposure rank',
       subtitle = paste0('Dataset ', example_seed, '; Spearman rho = ', round(predicted_rho, 3)),
       x = 'Actual standardized exposure rank', y = 'Predicted standardized exposure-rank score',
       caption = 'In-sample diagnostic; dashed line = agreement. Predictions are not pack-years.')
rank_limits <- range(c(example_patients$actual_exposure_rank_z, example_patients$predicted_exposure_rank_z))
p_predicted <- p_predicted + coord_equal(xlim = rank_limits, ylim = rank_limits)
save_plot(p_predicted, 'Predicted_vs_actual_smoking_rank')
p_raw <- ggplot(example_patients, aes(actual_exposure, predicted_exposure_rank_z)) +
  geom_point(alpha = 0.7, color = '#8361BB') +
  labs(title = 'Fitted exposure score vs recorded smoking',
       subtitle = paste0('Spearman rho = ', round(predicted_rho, 3)),
       x = paste0('Recorded exposure: ', exposure_column), y = 'Predicted standardized exposure-rank score',
       caption = 'Different units on the two axes; no raw smoking-amount prediction is made.')
save_plot(p_raw, 'Predicted_rank_vs_recorded_smoking')

# =============================================================================
# Supplement 1: model training losses
# =============================================================================
# The exporter does not save correlations at each epoch. Do not invent those curves.
training <- histories %>% filter(seed == example_seed) %>%
  pivot_longer(c(loss, exposure_loss, transition_loss), names_to = 'objective', values_to = 'value')
training_mean <- training %>% group_by(kind, objective, epoch) %>% summarise(value = mean(value), .groups = 'drop')
p_training <- ggplot(training, aes(epoch, value, color = kind)) +
  geom_line(aes(group = interaction(kind, repetition)), alpha = 0.18, linewidth = 0.3) +
  geom_line(data = training_mean, aes(group = kind), linewidth = 0.8) +
  facet_wrap(~ objective, scales = 'free_y') + scale_color_manual(values = kind_colors) +
  labs(title = 'Supplement 1: model optimization', x = 'Epoch', y = 'Training loss', color = 'Fit type',
       caption = 'Thin lines = individual fits; thick lines = means. No held-out validation curves.')
save_plot(p_training, 'Supp1_training_losses', 12, 4.5)

# =============================================================================
# Supplement 2: exposure-aware node priority and patient prevalence
# =============================================================================
if (has_truth) {
  p_node <- ggplot(example_nodes, aes(factor(set_membership, levels = c('A', 'C', 'B')), normalized_priority, fill = set_membership)) +
    geom_violin(alpha = 0.35, trim = TRUE) + geom_boxplot(width = 0.18, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.45) + scale_fill_manual(values = class_colors) +
    labs(title = 'Exposure-aware cohort node priority', x = 'Simulated node set', y = 'Normalized median-rank priority') + theme(legend.position = 'none')
  p_prevalence <- ggplot(example_nodes, aes(factor(set_membership, levels = c('A', 'C', 'B')), topk_patient_prevalence, fill = set_membership)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 0.8, alpha = 0.45) +
    scale_fill_manual(values = class_colors) + labs(title = 'Patient prevalence', x = 'Simulated node set', y = 'Fraction of patients selecting node') + theme(legend.position = 'none')
  save_plot(p_node + p_prevalence, 'Supp2_node_priority_and_prevalence', 11, 5)
} else {
  p_node <- example_nodes %>% arrange(median_patient_rank, node_idx) %>% slice_head(n = 30) %>%
    ggplot(aes(reorder(node_label, normalized_priority), normalized_priority)) + geom_col(fill = '#337A96') +
    coord_flip() + labs(title = 'Top cohort nodes', x = NULL, y = 'Normalized rank priority')
  save_plot(p_node, 'Top_cohort_nodes', 8, 8)
}

# =============================================================================
# Supplement 3: cancer-transition node score (paired embedding change)
# =============================================================================
if (has_truth) {
  p_cancer <- ggplot(example_nodes, aes(factor(set_membership, levels = c('A', 'C', 'B')), median_patient_transition_score, fill = set_membership)) +
    geom_violin(alpha = 0.35) + geom_boxplot(width = 0.18, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.45) + scale_fill_manual(values = class_colors) +
    labs(title = 'Supplement 3: paired-state node perturbation', x = 'Simulated node set',
         y = 'Median patient node-transition score',
         caption = 'Paired embedding-change magnitude; embeddings were learned with exposure supervision.') + theme(legend.position = 'none')
} else {
  p_cancer <- ggplot(example_nodes, aes(median_patient_transition_score)) + geom_histogram(bins = 40, fill = '#337A96', color = 'white') +
    labs(title = 'Paired-state node perturbation', x = 'Median patient node-transition score', y = 'Number of nodes')
}
save_plot(p_cancer, 'Supp3_node_transition_score')

# =============================================================================
# Supplement 4: selected edges and exposure-aware edge scores
# =============================================================================
edge_counts <- edges %>% group_by(seed, endpoint_group) %>% summarise(
  available = n(), selected_count = sum(selected), .groups = 'drop') %>% group_by(seed) %>%
  mutate(expected_count = sum(selected_count) * available / sum(available),
         selection_fraction = selected_count / available,
         fold_enrichment = ifelse(expected_count > 0, selected_count / expected_count, NA_real_)) %>% ungroup()
write_csv(edge_counts, file.path(output_dir, 'edge_selection_vs_graph_baseline.csv'))
p_edge_count <- edge_counts %>% filter(seed == example_seed) %>%
  ggplot(aes(endpoint_group, selected_count)) + geom_col(fill = '#ED7658', alpha = 0.8) +
  geom_point(aes(y = expected_count), shape = 95, size = 10) +
  labs(title = 'Selected directed edges by endpoint node classes', x = 'Endpoint node classes', y = 'Number selected',
       caption = 'Black marks = graph-composition expectation; endpoint classes are not edge ground truth.')
p_edge_score <- ggplot(example_edges, aes(endpoint_group, median_patient_priority, fill = endpoint_group)) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(title = 'Exposure-aware edge priority', x = 'Endpoint node classes', y = 'Median absolute score change on message suppression') + theme(legend.position = 'none')
save_plot(p_edge_count / p_edge_score, 'Supp4_edge_selection_and_priority', 10, 9)
p_edge_message <- ggplot(example_edges, aes(median_patient_message_transition_norm, median_patient_priority, color = endpoint_group)) +
  geom_point(alpha = 0.35, size = 1) + labs(title = 'Edge message change vs exposure priority',
    x = 'Median paired message-change norm', y = 'Median exposure edge priority', color = 'Endpoint classes')
save_plot(p_edge_message, 'Edge_message_change_vs_priority')

# =============================================================================
# Supplement 5: selected message-passing degree and weighted edge strength
# =============================================================================
# Directed incoming/outgoing counts stay separate. These are NOT biological rewiring.
degree_long <- degrees %>% filter(seed == example_seed) %>%
  pivot_longer(c(selected_incoming_degree, selected_outgoing_degree,
                 incoming_exposure_strength, outgoing_exposure_strength),
               names_to = 'measure', values_to = 'value')
p_degree <- ggplot(degree_long, aes(set_membership, value, fill = set_membership)) +
  geom_boxplot(outlier.alpha = 0.2) + facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  labs(title = 'Supplement 5: node connectivity in prioritized message channels',
       x = 'Node class', y = 'Patient-average count or strength',
       caption = 'Counts use patient top edges; strengths sum priorities over all incident channels.') + theme(legend.position = 'none')
if (has_truth) p_degree <- p_degree + scale_fill_manual(values = class_colors)
save_plot(p_degree, 'Supp5_selected_degree_and_edge_strength', 11, 7)
write_csv(degrees, file.path(output_dir, 'node_edge_summary_plot_data.csv'))

# =============================================================================
# Additional figure: patient-specific node priority heatmap
# =============================================================================
chosen <- example_nodes %>% arrange(median_patient_rank, node_idx) %>% slice_head(n = heatmap_top_nodes)
patient_order <- example_patients %>% arrange(actual_exposure, patient_id) %>% pull(patient_id)
heat <- patient_nodes %>% filter(seed == example_seed, node_idx %in% chosen$node_idx) %>%
  mutate(normalized_priority = (nrow(example_nodes) - median_rank) / max(nrow(example_nodes) - 1, 1),
         node_label = factor(node_label, levels = chosen$node_label),
         patient_id = factor(patient_id, levels = rev(patient_order)))
p_heat <- ggplot(heat, aes(node_label, patient_id, fill = normalized_priority)) + geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = 'Patient-specific contextual node priority', x = 'Top cohort nodes',
       y = 'Patients ordered by exposure', fill = 'Rank priority') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_plot(p_heat, 'Additional_patient_node_heatmap', 11, 7)

# =============================================================================
# Additional figures: patient-specific edges and attention changes
# =============================================================================
if (plot_patient_edges) {
  # Read only the needed columns to reduce memory use for large TCGA exports.
  patient_edges <- read_csv(file.path(datasets$analysis_dir[1], 'patient_edge_scores.csv.gz'),
    show_col_types = FALSE, progress = FALSE,
    col_types = cols_only(patient_id = col_character(), edge_id = col_integer(),
      median_rank = col_double(), median_priority = col_double(),
      median_attention_change = col_double(), median_message_transition_norm = col_double()))
  chosen_edges <- example_edges %>% arrange(median_patient_rank, edge_id) %>%
    slice_head(n = heatmap_top_edges) %>%
    mutate(edge_label = paste0(edge_id, ': ', source_label, ' -> ', target_label))
  edge_heat <- patient_edges %>% filter(edge_id %in% chosen_edges$edge_id) %>%
    left_join(chosen_edges %>% select(edge_id, edge_label), by = 'edge_id') %>%
    mutate(normalized_priority = (nrow(example_edges) - median_rank) / max(nrow(example_edges) - 1, 1),
      edge_label = factor(edge_label, levels = chosen_edges$edge_label),
      patient_id = factor(patient_id, levels = rev(patient_order)))
  p_edge_heat <- ggplot(edge_heat, aes(edge_label, patient_id, fill = normalized_priority)) +
    geom_tile() + scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = 'Patient-specific exposure edge priority', x = 'Directed edge',
      y = 'Patients ordered by exposure', fill = 'Rank priority') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  save_plot(p_edge_heat, 'Additional_patient_edge_heatmap', 11, 8)
  attention_summary <- patient_edges %>% group_by(edge_id) %>% summarise(
    median_patient_absolute_attention_change = median(abs(median_attention_change)), .groups = 'drop') %>%
    left_join(example_edges %>% select(edge_id, endpoint_group, median_patient_priority), by = 'edge_id')
  write_csv(attention_summary, file.path(output_dir, 'edge_attention_plot_data.csv'))
  p_attention <- ggplot(attention_summary, aes(median_patient_absolute_attention_change,
    median_patient_priority, color = endpoint_group)) + geom_point(alpha = 0.4, size = 1) +
    labs(title = 'Attention change vs exposure edge priority',
      x = 'Median patient absolute consensus attention change', y = 'Median exposure edge priority',
      color = 'Endpoint classes', caption = 'Attention change alone is not the exposure edge priority.')
  save_plot(p_attention, 'Additional_attention_change_vs_edge_priority')
  write_csv(edge_heat, file.path(output_dir, 'patient_edge_heatmap_data.csv'))
}

# =============================================================================
# Reviewer analysis: A versus C AUROC across independent simulated datasets
# =============================================================================
if (has_truth) {
  recovery <- nodes %>% group_by(seed) %>% summarise(
    `A vs C` = auc(set_membership, cohort_priority, 'A', 'C'),
    `A vs B` = auc(set_membership, cohort_priority, 'A', 'B'),
    `C vs B` = auc(set_membership, cohort_priority, 'C', 'B'), .groups = 'drop') %>%
    pivot_longer(-seed, names_to = 'comparison', values_to = 'AUROC')
  recovery$comparison <- factor(recovery$comparison, levels = c('A vs C', 'A vs B', 'C vs B'))
  recovery_summary <- recovery %>% group_by(comparison) %>% summarise(
    datasets = n(), mean = mean(AUROC), sd = sd(AUROC), min = min(AUROC), max = max(AUROC), .groups = 'drop')
  write_csv(recovery, file.path(output_dir, 'AUROC_by_dataset.csv'))
  write_csv(recovery_summary, file.path(output_dir, 'AUROC_summary.csv'))
  p_auc <- ggplot(recovery, aes(comparison, AUROC)) + geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey50') +
    geom_jitter(width = 0.08, height = 0, size = 2.5, color = '#337A96') +
    geom_errorbar(data = recovery_summary, aes(y = mean, ymin = mean - sd, ymax = mean + sd), width = 0.15, na.rm = TRUE) +
    geom_point(data = recovery_summary, aes(y = mean), shape = 95, size = 8) +
    coord_cartesian(ylim = c(0, 1)) + labs(title = 'Node recovery: A versus C is primary', x = NULL, y = 'AUROC',
      caption = 'Dots = independent datasets; bars = mean +/- SD when multiple datasets are supplied.')
  save_plot(p_auc, 'Reviewer_A_vs_C_recovery')

  # Top-node counts, expected counts and class-size-adjusted recall.
  top_nodes <- nodes %>% group_by(seed) %>% arrange(median_patient_rank, node_idx, .by_group = TRUE) %>%
    slice_head(n = cohort_top_k) %>% ungroup()
  totals <- nodes %>% count(seed, set_membership, name = 'set_size')
  top_counts <- top_nodes %>% count(seed, set_membership, name = 'selected')
  top_counts <- totals %>% left_join(top_counts, by = c('seed', 'set_membership')) %>%
    mutate(selected = replace_na(selected, 0L)) %>% group_by(seed) %>%
    mutate(k = sum(selected), expected = k * set_size / sum(set_size), recall = selected / set_size,
           fold_enrichment = selected / expected) %>% ungroup()
  write_csv(top_counts, file.path(output_dir, 'top_node_recovery.csv'))
  count_summary <- top_counts %>% group_by(set_membership) %>% summarise(
    mean_selected = mean(selected), expected = mean(expected), .groups = 'drop')
  p_top <- ggplot(count_summary, aes(set_membership, mean_selected, fill = set_membership)) +
    geom_col(alpha = 0.65) + geom_point(aes(y = expected), shape = 95, size = 10) +
    geom_point(data = top_counts, aes(y = selected), position = position_jitter(width = 0.08, height = 0)) +
    scale_fill_manual(values = class_colors) + labs(title = paste0('Top-', min(cohort_top_k, nrow(example_nodes)), ' cohort nodes'),
      x = 'Node set', y = 'Number selected', caption = 'Dots = independent datasets; black marks = random-selection expectation.') + theme(legend.position = 'none')
  save_plot(p_top, 'Reviewer_top_nodes_vs_expectation')
  p_recall <- ggplot(top_counts, aes(set_membership, recall, color = set_membership)) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.5) +
    scale_color_manual(values = class_colors) + coord_cartesian(ylim = c(0, 1)) +
    labs(title = 'Within-set recovery', x = 'Node set', y = 'Fraction recovered in cohort top list') + theme(legend.position = 'none')
  save_plot(p_recall, 'Reviewer_node_recall')

  # Simple text schematic, with counts from the actual illustrative dataset.
  n_a <- sum(example_nodes$set_membership == 'A'); n_c <- sum(example_nodes$set_membership == 'C')
  n_b <- sum(example_nodes$set_membership == 'B')
  schematic <- tibble(x = 1:3, y = 1, text = c(
    paste0(nrow(example_patients), ' paired patients\n', nrow(example_nodes), ' gene nodes\nRNA + methylation + CNV'),
    paste0('A: exposure-linked redistribution (', n_a, ')\nC: fixed transition shifts (', n_c, ')\nB: background (', n_b, ')'),
    'Attention and message passing\nPatient and cohort priorities\nInitialization / subsets / permutations'))
  p_design <- ggplot(schematic, aes(x, y, label = text)) + geom_label(size = 3.3, fill = '#F4F7FA') +
    xlim(0.4, 3.6) + ylim(0.7, 1.3) + theme_void() + labs(title = 'A. Paired multi-omics simulation')
  figure2 <- p_design / ((p_predicted + labs(title = 'B. Fitted exposure-rank score')) +
    (p_auc + labs(title = 'C. Node recovery')) + (p_top + labs(title = 'D. Top-ranked cohort nodes'))) +
    plot_layout(heights = c(0.5, 1))
  save_plot(figure2, 'Figure2_combined', 17, 9)
}

# =============================================================================
# Reviewer analysis: total change vs exposure and equal connectivity
# =============================================================================
audit_correlations <- audit %>% group_by(seed, measure) %>% summarise(
  exposure_rho = rho(magnitude, actual_exposure), .groups = 'drop')
write_csv(audit, file.path(output_dir, 'observed_magnitude_data.csv'))
write_csv(audit_correlations, file.path(output_dir, 'magnitude_exposure_correlations.csv'))
if (length(latent_list)) write_csv(bind_rows(latent_list), file.path(output_dir, 'latent_magnitude_ranges.csv'))
p_magnitude <- ggplot(audit_correlations, aes(measure, exposure_rho)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_jitter(width = 0.08, height = 0, size = 2) + coord_cartesian(ylim = c(-1, 1)) +
  labs(title = 'Total observed change vs exposure', x = NULL, y = 'Spearman correlation per dataset') +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
p_connectivity <- ggplot(input_degrees, aes(set_membership, degree)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 1) +
  labs(title = 'Actual input connectivity', x = 'Node set', y = 'Unique undirected degree')
save_plot(p_magnitude + p_connectivity, 'Reviewer_magnitude_and_connectivity', 12, 5)

# =============================================================================
# Confidence summaries: initialization, patient subsets and exposure permutation
# =============================================================================
# Frequencies and null fractions are descriptive. Do not combine into a p-value.
node_confidence <- example_nodes %>% select(node_label, set_membership,
  mean_initialization_stability, mean_subset_stability, mean_null_fraction) %>%
  pivot_longer(starts_with('mean_'), names_to = 'measure', values_to = 'value')
p_confidence <- ggplot(node_confidence, aes(set_membership, value, fill = set_membership)) +
  geom_boxplot(outlier.alpha = 0.2) + facet_wrap(~ measure) +
  labs(title = 'Node-priority confidence summaries', x = 'Node set', y = 'Frequency / descriptive null fraction',
       caption = 'Higher stability and lower null fraction have different interpretations.') + theme(legend.position = 'none')
if (has_truth) p_confidence <- p_confidence + scale_fill_manual(values = class_colors)
save_plot(p_confidence, 'Additional_node_confidence', 12, 4.5)
edge_confidence <- example_edges %>% select(endpoint_group, mean_initialization_stability,
  mean_subset_stability, mean_null_fraction) %>%
  pivot_longer(starts_with('mean_'), names_to = 'measure', values_to = 'value')
p_edge_confidence <- ggplot(edge_confidence, aes(endpoint_group, value)) + geom_boxplot(outlier.alpha = 0.1) +
  facet_wrap(~ measure) + labs(title = 'Edge-priority confidence summaries', x = 'Endpoint node classes',
    y = 'Frequency / descriptive null fraction', caption = 'Endpoint classes describe graph composition, not edge-recovery truth.')
save_plot(p_edge_confidence, 'Additional_edge_confidence', 12, 4.5)

# Final per-run fitted correlation; not an epoch-by-epoch validation curve.
run_correlations <- run_patients %>% group_by(seed, kind, `repeat`) %>% summarise(
  fitted_exposure_rho = rho(predicted_exposure_rank_score[included_in_fit], supplied_exposure[included_in_fit]),
  .groups = 'drop')
write_csv(run_correlations, file.path(output_dir, 'final_fit_correlations.csv'))
p_fit_corr <- run_correlations %>% filter(seed == example_seed) %>%
  ggplot(aes(kind, fitted_exposure_rho, color = kind)) + geom_jitter(width = 0.1, height = 0, size = 2) +
  scale_color_manual(values = kind_colors) + labs(title = 'Final fitted exposure alignment across repetitions',
    x = 'Fit type', y = 'In-sample Spearman correlation',
    caption = 'Permutation runs use their shuffled targets. High fitted correlation is not evidence of specificity.') + theme(legend.position = 'none')
save_plot(p_fit_corr, 'Additional_final_fit_correlations')

# =============================================================================
# Record scope and reproducibility
# =============================================================================
writeLines(c(
  'All numerical results were read from the supplied extension outputs.',
  paste('Illustrative dataset:', example_seed),
  paste('Number of datasets:', nrow(datasets)),
  'Single-dataset distributions do not establish uncertainty across simulations.',
  'Patient predictions are standardized exposure ranks, not raw smoking amounts.',
  'No p-values or node-wise significance tests were computed.',
  'Edge classes refer to endpoint node sets, not independently validated edge truth.',
  'Attention/message changes on fixed topology are not biological rewiring.',
  'No CpG-to-gene relationship validation is performed.',
  'Per-epoch correlation and held-out validation are unavailable in these exports.',
  'Latent magnitude checks, when available, are read from simulator diagnostics.'
), file.path(output_dir, 'analysis_notes.txt'))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'R_session_info.txt'))
message('Saved plots and supporting tables to: ', normalizePath(output_dir))

# #############################################################################
# PART 2: TCGA-specific plots (formerly plot_tcga.R)
# #############################################################################
# =============================================================================
# Settings for TCGA-specific descriptive panels
# =============================================================================
node_fraction <- 0.10
edge_fraction <- fromJSON(file.path(analysis_dir, 'extended_score_metadata.json'))$edge_top_fraction
patient_list <- patients$patient_id

# =============================================================================
# Recorded exposure groups: keep tied exposures together
# =============================================================================
# Recorded zero is not relabelled as never-smoking. Thresholds are descriptive.
cutpoints <- quantile(patients$actual_exposure, probs = c(1/3, 2/3), names = FALSE)
if (cutpoints[1] == cutpoints[2]) {
  message('Tied tertile thresholds: using two recorded-exposure groups instead of splitting ties.')
  patients <- patients %>% mutate(exposure_group = if_else(actual_exposure <= cutpoints[1], 'Lower recorded exposure', 'Higher recorded exposure'))
} else {
  patients <- patients %>% mutate(exposure_group = case_when(
    actual_exposure <= cutpoints[1] ~ 'Low',
    actual_exposure <= cutpoints[2] ~ 'Medium', TRUE ~ 'High'))
}
patients$exposure_group <- factor(patients$exposure_group,
  levels = if (cutpoints[1] == cutpoints[2]) c('Lower recorded exposure','Higher recorded exposure') else c('Low','Medium','High'))
write_csv(patients, file.path(output_dir, 'TCGA_patient_plot_data.csv'))
groups <- patients %>% group_by(exposure_group) %>% summarise(
  n_patients = n(), min_exposure = min(actual_exposure), max_exposure = max(actual_exposure),
  median_exposure = median(actual_exposure), .groups = 'drop')
write_csv(groups, file.path(output_dir, 'TCGA_exposure_group_definitions.csv'))

# =============================================================================
# Available cohort characteristics and missingness in supplied metadata
# =============================================================================
clinical <- read_table(file.path(data_dir, 'metadata_samples.csv')) %>%
  mutate(patient_id = as.character(patient_id)) %>% distinct(patient_id, .keep_all = TRUE) %>%
  filter(patient_id %in% patients$patient_id)
write_csv(clinical, file.path(output_dir, 'TCGA_supplied_patient_metadata.csv'))
cohort_counts <- tibble(characteristic = 'Analyzed patients', value = 'All', count = nrow(clinical))
for (field in intersect(c('smoking_status','gender','histology','batch'), names(clinical))) {
  tab <- clinical %>% transmute(value = replace_na(as.character(.data[[field]]), 'Missing in export')) %>%
    count(value, name = 'count') %>% mutate(characteristic = field) %>% select(characteristic,value,count)
  cohort_counts <- bind_rows(cohort_counts, tab)
}
write_csv(cohort_counts, file.path(output_dir, 'TCGA_available_cohort_counts.csv'))
cohort_numeric <- tibble()
for (field in intersect(c('age','smoking_intensity','tumor_purity'), names(clinical))) {
  x <- suppressWarnings(as.numeric(clinical[[field]])); valid <- is.finite(x)
  cohort_numeric <- bind_rows(cohort_numeric, tibble(variable = field, n_observed_in_export = sum(valid),
    n_missing_in_export = sum(!valid), median = if (any(valid)) median(x[valid]) else NA_real_,
    q25 = if (any(valid)) unname(quantile(x[valid], .25)) else NA_real_,
    q75 = if (any(valid)) unname(quantile(x[valid], .75)) else NA_real_))
}
write_csv(cohort_numeric, file.path(output_dir, 'TCGA_available_cohort_numeric.csv'))

# =============================================================================
# Figure 3B: patient-specific node contribution among prioritized nodes
# =============================================================================
# Raw contribution magnitude is used; median rank of a fixed top list is constrained.
node_k <- max(1, ceiling(node_fraction * nrow(nodes)))
patient_selected_nodes <- patient_nodes %>% group_by(patient_id) %>%
  arrange(median_rank, node_idx, .by_group = TRUE) %>% slice_head(n = node_k) %>% ungroup()
node_summary <- patient_selected_nodes %>% group_by(patient_id) %>%
  summarise(top_node_median_priority = median(median_priority),
            top_node_mean_priority = mean(median_priority), .groups = 'drop') %>%
  left_join(patients %>% select(patient_id, exposure_group, actual_exposure), by = 'patient_id')
p_group_nodes <- ggplot(node_summary, aes(exposure_group, top_node_median_priority, fill = exposure_group)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .12, height = 0, size = 1.5, alpha = .7) +
  labs(title = 'B. Patient-specific prioritized node contributions',
       subtitle = paste0('Top ', round(node_fraction*100), '% of nodes per patient (n = ',node_k,')'),
       x = 'Recorded exposure group', y = 'Median absolute node contribution',
       caption = 'Descriptive fitted-score distributions; no group-test p-values.') + theme(legend.position = 'none')
save_plot(p_group_nodes, 'Figure3B_patient_node_contributions', 8, 5)

# =============================================================================
# Figure 3C: incoming edge-priority strength among prioritized nodes
# =============================================================================
node_edges <- read_table(file.path(analysis_dir, 'patient_node_edge_summary.csv.gz')) %>%
  mutate(patient_id = as.character(patient_id))
selected_node_edges <- patient_selected_nodes %>% select(patient_id,node_idx) %>%
  left_join(node_edges, by = c('patient_id','node_idx'))
network_summary <- selected_node_edges %>% group_by(patient_id) %>% summarise(
  prioritized_node_incoming_strength = mean(incoming_exposure_strength),
  prioritized_node_selected_incoming_degree = mean(selected_incoming_degree),
  prioritized_node_selected_outgoing_degree = mean(selected_outgoing_degree), .groups = 'drop') %>%
  left_join(patients %>% select(patient_id, exposure_group, actual_exposure), by = 'patient_id')
p_group_edges <- ggplot(network_summary, aes(exposure_group, prioritized_node_incoming_strength, fill = exposure_group)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .12, height = 0, size = 1.5, alpha = .7) +
  labs(title = 'C. Edge-priority strength at prioritized nodes', x = 'Recorded exposure group',
       y = 'Mean incoming exposure-edge strength',
       caption = 'Strength sums all incoming edge priorities at each selected node; not biological rewiring.') + theme(legend.position = 'none')
save_plot(p_group_edges, 'Figure3C_patient_edge_strength', 8, 5)
patient_summaries <- node_summary %>% left_join(network_summary %>%
  select(-exposure_group,-actual_exposure), by = 'patient_id')
write_csv(patient_summaries, file.path(output_dir, 'TCGA_patient_node_edge_summaries.csv'))
# Figure 3A is transition magnitude, NOT the fitted exposure prediction.
figure3abc <- (p_transition + labs(title = 'A. Paired transition magnitude vs exposure')) /
  (p_group_nodes + p_group_edges) + plot_layout(heights = c(1, 1))
save_plot(figure3abc, 'Figure3_ABC_no_pathway_panel', 14, 10)

# =============================================================================
# Supplement 7: cohort node and edge score distributions and selection boundaries
# =============================================================================
# The boundary is a descriptive top-fraction ranking cutoff, not significance.
# Use genes for the gene panel; individual CpG nodes have separate panels below.
genes <- nodes %>% filter(grepl('^GENE:', node_label)) %>% arrange(median_patient_rank,node_idx)
gene_k <- max(1, ceiling(node_fraction*nrow(genes)))
gene_threshold <- genes$normalized_priority[gene_k]
edge_order <- edges %>% arrange(median_patient_rank,edge_id)
edge_k <- max(1, ceiling(edge_fraction*nrow(edge_order)))
edge_threshold <- edge_order$normalized_priority[edge_k]
p_node_dist <- ggplot(genes, aes(normalized_priority)) + geom_histogram(bins = 40, fill = '#337A96', color = 'white') +
  geom_vline(xintercept = gene_threshold, linetype = 'dashed') +
  labs(title = 'Gene-node priority distribution', x = 'Normalized cohort rank priority', y = 'Gene count')
p_edge_dist <- ggplot(edges, aes(normalized_priority)) + geom_histogram(bins = 40, fill = '#ED7658', color = 'white') +
  geom_vline(xintercept = edge_threshold, linetype = 'dashed') +
  labs(title = 'Directed-edge priority distribution', x = 'Normalized cohort rank priority', y = 'Directed edge count')
save_plot(p_node_dist + p_edge_dist, 'TCGA_Supp7_priority_distributions', 12, 5)
write_csv(tibble(kind = c('gene','directed_edge'), fraction = c(node_fraction,edge_fraction),
  k = c(gene_k,edge_k), normalized_priority_boundary = c(gene_threshold,edge_threshold),
  tie_rule = 'Sort by median patient rank, then node_idx or edge_id; boundary ties may extend beyond k'),
  file.path(output_dir, 'TCGA_selection_boundaries.csv'))
write_csv(genes %>% slice_head(n = gene_k), file.path(output_dir, 'TCGA_top10percent_genes.csv'))
write_csv(edge_order %>% slice_head(n = edge_k), file.path(output_dir, 'TCGA_top_edges.csv'))

# =============================================================================
# Supplement 8: node contribution vs network strength, per patient
# =============================================================================
node_edge_rho <- rho(patient_summaries$top_node_median_priority,
                    patient_summaries$prioritized_node_incoming_strength)
p_relationship <- ggplot(patient_summaries, aes(prioritized_node_incoming_strength,
  top_node_median_priority, color = exposure_group)) + geom_point(size = 2, alpha = .8) +
  labs(title = 'TCGA Supplement 8: patient node and edge summaries',
    subtitle = paste0('Descriptive Spearman rho = ', round(node_edge_rho,3)),
    x = 'Mean incoming edge-priority strength at prioritized nodes',
    y = 'Median prioritized-node contribution', color = 'Recorded exposure group',
    caption = 'Both summaries come from the same fitted model; correlation is not independent validation.')
save_plot(p_relationship, 'TCGA_Supp8_node_edge_relationship', 8, 5)
write_csv(tibble(comparison='patient_node_vs_edge_strength',spearman_rho=node_edge_rho),
  file.path(output_dir,'TCGA_node_edge_correlation.csv'))

# =============================================================================
# Supplement 10: CpG priorities and candidate CpG-gene message channels
# =============================================================================
cpg_nodes <- nodes %>% filter(grepl('^CpG:',node_label)) %>% arrange(median_patient_rank,node_idx)
cpg_edges <- edges %>% filter(edge_type %in% c('cpg_to_gene','gene_to_cpg')) %>% arrange(median_patient_rank,edge_id)
if (nrow(cpg_nodes)>0 && nrow(cpg_edges)>0) {
  cpg_k <- max(1,ceiling(.10*nrow(cpg_nodes)))
  cpg_edge_k <- max(1,ceiling(.10*nrow(cpg_edges)))
  p_cpg <- ggplot(cpg_nodes,aes(normalized_priority)) + geom_histogram(bins=25,fill='#8361BB',color='white') +
    geom_vline(xintercept=cpg_nodes$normalized_priority[cpg_k],linetype='dashed') +
    labs(title='CpG node priority',x='Normalized cohort rank priority',y='CpG count')
  p_cpg_edge <- ggplot(cpg_edges,aes(normalized_priority)) + geom_histogram(bins=30,fill='#337A96',color='white') +
    geom_vline(xintercept=cpg_edges$normalized_priority[cpg_edge_k],linetype='dashed') +
    labs(title='Candidate CpG-gene channel priority',x='Normalized cohort rank priority',y='Directed channel count')
  selected_cpgs <- cpg_nodes %>% slice_head(n=cpg_k)
  cpg_heat <- patient_nodes %>% filter(node_idx %in% selected_cpgs$node_idx) %>%
    mutate(node_label=factor(node_label,levels=rev(selected_cpgs$node_label)),
      patient_id=factor(patient_id,levels=patient_order),
      rank_priority=(nrow(nodes)-median_rank)/max(1,nrow(nodes)-1))
  p_cpg_heat <- ggplot(cpg_heat,aes(patient_id,node_label,fill=rank_priority)) + geom_tile() +
    scale_fill_viridis_c(limits=c(0,1)) + labs(title='Patient-specific priorities of top CpG nodes',
      x='Patients ordered by recorded exposure',y='CpG identifier',fill='Rank priority',
      caption='Prior CpG-gene connections do not establish regulation; scores do not indicate methylation direction.') +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=7))
  save_plot((p_cpg_edge+p_cpg)/p_cpg_heat,'TCGA_Supp10_CpG_summaries',13,9)
  write_csv(selected_cpgs,file.path(output_dir,'TCGA_top_CpG_nodes.csv'))
  write_csv(cpg_edges %>% slice_head(n=cpg_edge_k),file.path(output_dir,'TCGA_top_CpG_gene_edges.csv'))
} else {
  message('No CpG-node/CpG-gene exports: skipping Supplement 10.')
}


# =============================================================================
# Separate top gene and CpG nodes
# =============================================================================

for (node_type in c("GENE", "CpG")) {

  # Keep the original scores; only filter which node type is displayed.
  ranked_nodes <- example_nodes %>%
    filter(startsWith(node_label, paste0(node_type, ":"))) %>%
    arrange(median_patient_rank, node_idx)

  if (nrow(ranked_nodes) == 0) next

  top_nodes <- ranked_nodes %>%
    slice_head(n = 30) %>%
    mutate(node_label = factor(node_label, levels = rev(node_label)))

  plot_title <- if (node_type == "GENE") {
    "Top cohort gene nodes"
  } else {
    "Top cohort CpG nodes"
  }

  plot_color <- if (node_type == "GENE") "#337A96" else "#8361BB"

  p_top <- ggplot(top_nodes, aes(node_label, normalized_priority)) +
    geom_col(fill = plot_color, width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Weighted node score",
      caption = "Original cohort scores; top 30 nodes within the indicated node type."
    ) +
    theme(plot.title = element_text(face = "bold"))

  output_name <- if (node_type == "GENE") {
    "Top_cohort_gene_nodes"
  } else {
    "Top_cohort_CpG_nodes"
  }

  save_plot(p_top, output_name, 8, 8)

  write_csv(
    ranked_nodes,
    file.path(output_dir, paste0(output_name, "_all_ranked.csv"))
  )
}

# =============================================================================
# Supplementary Figures 6–10
# Paste at the END of plot_tcga.R
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

set.seed(42)

fig_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

group_colors <- c(
  Low = "#55964A",
  Medium = "#9773AF",
  High = "#FF4040",
  "Lower recorded exposure" = "#55964A",
  "Higher recorded exposure" = "#FF4040"
)

patient_info <- patients %>%
  select(patient_id, actual_exposure, exposure_group) %>%
  arrange(exposure_group, actual_exposure, patient_id)


# =============================================================================
# Supplement 6: model optimization
# Thin lines = individual fits; thick lines = means
# No held-out validation curves are available
# =============================================================================

h6 <- histories %>%
  filter(seed == example_seed) %>%
  pivot_longer(
    c(exposure_loss, transition_loss),
    names_to = "objective",
    values_to = "value"
  ) %>%
  mutate(
    objective = factor(
      objective,
      levels = c("exposure_loss", "transition_loss"),
      labels = c(
        "Exposure-rank loss",
        "Transition-reconstruction loss"
      )
    )
  )

h6_mean <- h6 %>%
  group_by(kind, objective, epoch) %>%
  summarise(value = mean(value), .groups = "drop")

p6 <- ggplot(h6, aes(epoch, value, color = kind)) +
  geom_line(
    aes(group = interaction(kind, repetition)),
    alpha = 0.15,
    linewidth = 0.3
  ) +
  geom_line(
    data = h6_mean,
    aes(group = kind),
    linewidth = 0.9
  ) +
  facet_wrap(~ objective, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c(
      initialization = "#2886B5",
      subset = "#C17A26",
      permutation = "#A14792"
    )
  ) +
  labs(
    title = "Model optimization in TCGA",
    x = "Epoch",
    y = "Training loss",
    color = "Fit type",
    caption = paste(
      "Thin lines: individual fits; bold lines: means.",
      "Training diagnostics, not validation."
    )
  ) +
  fig_theme

save_plot(p6, "TCGA_Supp6_training", 8, 8)


# =============================================================================
# Helper: score distribution with histogram, density and top-10% cutoff
# =============================================================================

score_distribution <- function(
    values, title, x_label, line_color = "red") {

  values <- values[is.finite(values)]
  if (!length(values)) stop("No finite scores for ", title)

  cutoff <- sort(values, decreasing = TRUE)[
    max(1, ceiling(0.10 * length(values)))
  ]

  p <- ggplot(tibble(score = values), aes(score)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      fill = "grey75",
      color = "grey30"
    ) +
    geom_vline(
      xintercept = cutoff,
      linetype = "dashed",
      color = line_color
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("Top 10% cutoff = ", signif(cutoff, 3)),
      hjust = 1.05,
      vjust = 1.5,
      size = 3.5,
      color = line_color
    ) +
    labs(title = title, x = x_label, y = "Density") +
    fig_theme

  if (length(unique(values)) > 1) {
    p <- p + geom_density(
      color = "black",
      linewidth = 0.7,
      trim = TRUE
    )
  }

  p
}


# =============================================================================
# Supplement 7: cohort gene-node and edge score distributions
# Original normalized cohort scores are retained
# =============================================================================

gene7 <- nodes %>%
  filter(startsWith(node_label, "GENE:")) %>%
  arrange(median_patient_rank, node_idx)

edge7 <- edges %>%
  arrange(median_patient_rank, edge_id)

p7a <- score_distribution(
  gene7$normalized_priority,
  "Distribution of weighted gene-node scores",
  "Weighted node score"
)

p7b <- score_distribution(
  edge7$normalized_priority,
  "Distribution of weighted edge scores",
  "Weighted edge score"
)

save_plot(
  (p7a / p7b) + plot_annotation(tag_levels = "A"),
  "TCGA_Supp7_distributions_revised",
  8, 9
)

write_csv(
  tibble(
    type = c("Gene", "Directed edge"),
    selected_n = c(
      ceiling(0.10 * nrow(gene7)),
      ceiling(0.10 * nrow(edge7))
    ),
    cutoff = c(
      gene7$normalized_priority[ceiling(0.10 * nrow(gene7))],
      edge7$normalized_priority[ceiling(0.10 * nrow(edge7))]
    )
  ),
  file.path(output_dir, "TCGA_Supp7_cutoffs.csv")
)


# =============================================================================
# Supplement 8: selected-edge strength versus patient-specific node score
#
# 1. Select top 10% edges separately for each patient.
# 2. Sum their raw edge priorities at each target node.
# 3. Select top 10% nodes by this incoming strength.
# 4. Compare mean strength with median node score in those SAME nodes.
#
# Directed edges are retained. Ties are broken using edge/node identifiers.
# =============================================================================

if (!exists("patient_edges")) {
  patient_edges <- read_csv(
    file.path(analysis_dir, "patient_edge_scores.csv.gz"),
    show_col_types = FALSE,
    col_types = cols_only(
      patient_id = col_character(),
      edge_id = col_integer(),
      median_rank = col_double(),
      median_priority = col_double()
    )
  )
}

edge_k8 <- ceiling(0.10 * nrow(edges))
node_k8 <- ceiling(0.10 * nrow(nodes))

selected8 <- patient_edges %>%
  group_by(patient_id) %>%
  arrange(median_rank, edge_id, .by_group = TRUE) %>%
  slice_head(n = edge_k8) %>%
  ungroup() %>%
  left_join(
    edges %>% select(edge_id, target_idx),
    by = "edge_id"
  )

strength8 <- selected8 %>%
  group_by(patient_id, target_idx) %>%
  summarise(
    selected_incoming_strength = sum(median_priority),
    .groups = "drop"
  ) %>%
  rename(node_idx = target_idx)

node8 <- patient_nodes %>%
  select(patient_id, node_idx, median_priority) %>%
  left_join(strength8, by = c("patient_id", "node_idx")) %>%
  mutate(
    selected_incoming_strength =
      replace_na(selected_incoming_strength, 0)
  ) %>%
  group_by(patient_id) %>%
  arrange(
    desc(selected_incoming_strength),
    node_idx,
    .by_group = TRUE
  ) %>%
  slice_head(n = node_k8) %>%
  summarise(
    mean_strength = mean(selected_incoming_strength),
    median_node_score = median(median_priority),
    .groups = "drop"
  ) %>%
  left_join(patient_info, by = "patient_id")

r8 <- if (
  sd(node8$mean_strength) > 0 &&
  sd(node8$median_node_score) > 0
) {
  cor(
    node8$mean_strength,
    node8$median_node_score,
    method = "spearman"
  )
} else {
  NA_real_
}

p8 <- ggplot(node8, aes(mean_strength, median_node_score)) +
  geom_point(
    aes(color = exposure_group),
    size = 2,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "#4278FF",
    fill = "grey65",
    se = TRUE
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Selected-edge strength vs node score",
    subtitle = paste0(
      "Descriptive Spearman rho = ", round(r8, 3)
    ),
    x = "Mean selected incoming edge strength",
    y = "Median patient-specific node score",
    color = "Recorded exposure group",
    caption = paste(
      "Top 10% edges per patient; top 10% nodes by incoming strength.",
      "Both scores are model-derived."
    )
  ) +
  fig_theme

save_plot(p8, "TCGA_Supp8_selected_edge_strength", 8, 6)

write_csv(
  node8,
  file.path(output_dir, "TCGA_Supp8_patient_values.csv")
)


# =============================================================================
# Helper: patient heatmap
# Patients ordered by recorded exposure within groups
# Raw node contributions or within-node z-scores are displayed
# =============================================================================

patient_heatmap <- function(selection, title, row_scale = FALSE) {

  h <- patient_nodes %>%
    inner_join(selection, by = "node_idx") %>%
    left_join(patient_info, by = "patient_id") %>%
    group_by(node_idx) %>%
    mutate(
      display_score = if (row_scale) {
        if (sd(median_priority) > 0) {
          (median_priority - mean(median_priority)) /
            sd(median_priority)
        } else {
          rep(0, n())
        }
      } else {
        median_priority
      }
    ) %>%
    ungroup() %>%
    mutate(
      patient_id = factor(
        patient_id,
        levels = patient_info$patient_id
      ),
      display_label = factor(
        display_label,
        levels = rev(selection$display_label)
      )
    )

  p <- ggplot(
    h,
    aes(patient_id, display_label, fill = display_score)
  ) +
    geom_tile() +
    facet_grid(
      ~ exposure_group,
      scales = "free_x",
      space = "free_x"
    ) +
    labs(
      title = title,
      x = "Patients ordered by recorded exposure",
      y = NULL
    ) +
    fig_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )

  if (row_scale) {
    p <- p + scale_fill_gradient2(
      low = "#25259A",
      mid = "white",
      high = "#E52727",
      midpoint = 0,
      name = "Within-node\nz-score"
    )
  } else {
    p <- p + scale_fill_gradientn(
      colors = c("#76E600", "#FFE000", "#FF3B00"),
      name = "Patient-specific\nnode score"
    )
  }

  list(plot = p, data = h)
}


# =============================================================================
# Supplement 9: selected smoking/lung-cancer pathway genes
#
# Provide pathway_genes.csv in GNN_Input_Data with these columns:
# node_label,gene_symbol
#
# node_label must match the model's GENE:<Entrez ID>.
# Use your verified pathway intersection or union.
# The figure is skipped if the file is absent.
# =============================================================================

pathway_file <- file.path(data_dir, "pathway_genes.csv")

if (file.exists(pathway_file)) {

  pathway9 <- read_csv(
    pathway_file,
    show_col_types = FALSE
  )

  stopifnot(
    all(c("node_label", "gene_symbol") %in% names(pathway9))
  )

  pathway9 <- pathway9 %>%
    distinct(node_label, .keep_all = TRUE)

  missing9 <- anti_join(pathway9, gene7, by = "node_label")

  write_csv(
    missing9,
    file.path(output_dir, "TCGA_Supp9_unmatched_genes.csv")
  )

  selection9 <- pathway9 %>%
    inner_join(
      gene7 %>% select(node_label, node_idx),
      by = "node_label"
    ) %>%
    transmute(
      node_idx,
      display_label = make.unique(as.character(gene_symbol))
    )

  if (nrow(selection9) > 0) {

    h9 <- patient_heatmap(
      selection9,
      "Patient-specific scores of selected pathway genes"
    )

    save_plot(
      h9$plot,
      "TCGA_Supp9_pathway_gene_heatmap",
      11,
      max(6, 0.18 * nrow(selection9) + 2)
    )

    write_csv(
      h9$data,
      file.path(output_dir, "TCGA_Supp9_heatmap_values.csv")
    )
  }

} else {
  message(
    "Supplement 9 skipped: supply pathway_genes.csv ",
    "with node_label and gene_symbol."
  )
}


# =============================================================================
# Supplement 10: candidate CpG-gene edges and CpG node priorities
# A: all candidate CpG-gene directed-edge scores
# B: all CpG cohort node scores
# C: top 10% CpGs, with patient-specific scores standardized within each CpG
# =============================================================================

cpg10 <- nodes %>%
  filter(startsWith(node_label, "CpG:")) %>%
  arrange(median_patient_rank, node_idx)

reg10 <- edges %>%
  filter(edge_type %in% c("cpg_to_gene", "gene_to_cpg")) %>%
  arrange(median_patient_rank, edge_id)

if (nrow(cpg10) > 0 && nrow(reg10) > 0) {

  p10a <- score_distribution(
    reg10$normalized_priority,
    "Candidate CpG-gene edge scores",
    "Weighted edge score",
    "black"
  )

  p10b <- score_distribution(
    cpg10$normalized_priority,
    "CpG weighted node scores",
    "Weighted node score",
    "black"
  )

  selected10 <- cpg10 %>%
    slice_head(n = ceiling(0.10 * nrow(cpg10)))

  selection10 <- selected10 %>%
    transmute(
      node_idx,
      display_label = sub("^CpG:", "", node_label)
    )

  h10 <- patient_heatmap(
    selection10,
    "Patient-specific scores of top 10% CpG nodes",
    row_scale = TRUE
  )

  save_plot(
    ((p10a + p10b) / h10$plot) +
      plot_layout(heights = c(1, 1.4)) +
      plot_annotation(tag_levels = "A"),
    "TCGA_Supp10_CpG_revised",
    13, 10
  )

  write_csv(
    h10$data,
    file.path(output_dir, "TCGA_Supp10_heatmap_values.csv")
  )

  write_csv(
    selected10,
    file.path(output_dir, "TCGA_Supp10_selected_CpGs.csv")
  )
}


# #############################################################################
# NEW SECTIONS
# #############################################################################
run_config <- fromJSON(file.path(analysis_dir, 'config.json'))
run_label <- paste0(nrow(example_patients), ' patients, ', run_config$epochs, ' epochs')

# =============================================================================
# N1. Fit check: in-sample fit shown next to the shuffled-label fits
# =============================================================================
# The left panel shows the exposure objective trained; the right panel shows that
# fit alone is not evidence of exposure specificity (shuffled labels also fit).
fit_runs <- run_correlations %>% filter(seed == example_seed,
                                         kind %in% c('initialization', 'permutation')) %>%
  mutate(labels = factor(if_else(kind == 'initialization', 'Real labels\n(initializations)',
                                 'Shuffled labels\n(permutations)'),
                         levels = c('Real labels\n(initializations)', 'Shuffled labels\n(permutations)')))
fit_summary <- fit_runs %>% group_by(labels) %>%
  summarise(n_fits = n(), median_rho = median(fitted_exposure_rho),
            min_rho = min(fitted_exposure_rho), max_rho = max(fitted_exposure_rho), .groups = 'drop')
write_csv(fit_summary %>% mutate(labels = gsub('\\n', ' ', labels)), file.path(output_dir, 'Fit_check_summary.csv'))
p_fit_left <- p_predicted + labs(title = 'A. Exposure objective was learned',
  subtitle = paste0(run_label, '; in-sample Spearman rho = ', round(predicted_rho, 2)),
  caption = 'Median fitted score across initializations. In-sample fit, not prediction.')
p_fit_right <- ggplot(fit_runs, aes(labels, fitted_exposure_rho, color = labels)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, color = 'grey40') +
  geom_jitter(width = 0.08, height = 0, size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c('#2886B5', '#A14792')) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = 'B. Fit with real vs shuffled labels', x = NULL,
       y = 'In-sample Spearman (fitted vs supplied)',
       caption = 'The model fits any supplied ranking; evidence for exposure nodes\ncomes from node-level stability and real-vs-shuffled differences.') +
  theme(legend.position = 'none')
save_plot(p_fit_left + p_fit_right + plot_layout(widths = c(1, 1)), 'Fit_check_real_vs_shuffled', 13, 6)

# =============================================================================
# N2. Real-vs-shuffled node agreement
# =============================================================================
# Uses each run's saved patient-by-node |contribution| (attributions.npz 'score').
# Per patient: rank agreement between pairs of runs.
#   real-real     = initialization vs initialization (seed-to-seed variability)
#   real-shuffled = initialization vs permutation
# If real-real is clearly higher, highlighted nodes depend on the real labels and
# not only on the random start. Descriptive; no p-values.
read_npy_from_npz <- function(npz, key) {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(npz, files = paste0(key, '.npy'), exdir = tmp)
  path <- file.path(tmp, paste0(key, '.npy'))
  bytes <- readBin(path, 'raw', file.info(path)$size)
  if (!identical(bytes[1:6], as.raw(c(0x93, 0x4E, 0x55, 0x4D, 0x50, 0x59)))) stop('Not a .npy file: ', key)
  major <- as.integer(bytes[7])
  if (major == 1) { hlen <- sum(as.integer(bytes[9:10]) * 256^(0:1)); start <- 11 }
  else { hlen <- sum(as.integer(bytes[9:12]) * 256^(0:3)); start <- 13 }
  header <- rawToChar(bytes[start:(start + hlen - 1)])
  descr <- sub(".*'descr': *'([^']+)'.*", '\\1', header)
  if (grepl("'fortran_order': *True", header)) stop('Fortran-ordered arrays are not supported')
  shape <- as.integer(strsplit(gsub('[^0-9,]', '', sub(".*'shape': *\\(([^)]*)\\).*", '\\1', header)), ',')[[1]])
  shape <- shape[!is.na(shape)]
  size <- switch(descr, '<f4' = 4, '<f8' = 8, stop('Unsupported dtype ', descr))
  values <- readBin(bytes[(start + hlen):length(bytes)], 'numeric', n = prod(shape), size = size, endian = 'little')
  if (length(shape) == 2) matrix(values, nrow = shape[1], ncol = shape[2], byrow = TRUE) else values
}

agreement_runs <- read_table(file.path(analysis_dir, 'runs.csv')) %>%
  filter(kind %in% c('initialization', 'permutation')) %>%
  mutate(folder = file.path(analysis_dir, paste0(kind, '_', .data[['repeat']])))
has_npz <- nrow(agreement_runs) > 0 && all(file.exists(file.path(agreement_runs$folder, 'attributions.npz')))
if (has_npz && sum(agreement_runs$kind == 'initialization') >= 2 &&
    any(agreement_runs$kind == 'permutation')) {
  score_list <- lapply(agreement_runs$folder, function(f) read_npy_from_npz(file.path(f, 'attributions.npz'), 'score'))
  run_patient_ids <- as.character(read_table(file.path(agreement_runs$folder[1], 'fit_membership.csv'))$patient_id)
  n_nodes <- ncol(score_list[[1]])
  top_k <- min(run_config$top_k, n_nodes)
  is_init <- agreement_runs$kind == 'initialization'; is_perm <- !is_init
  pair_mean <- function(m, a, b, same) {
    block <- m[a, b, drop = FALSE]
    if (same) block[upper.tri(block)] else as.vector(block)
  }
  patient_agreement <- list()
  cohort_rank_by_run <- matrix(NA_real_, nrow(agreement_runs), n_nodes)
  rank_sum <- vector('list', nrow(agreement_runs))
  for (r in seq_len(nrow(agreement_runs))) rank_sum[[r]] <- matrix(NA_real_, length(run_patient_ids), n_nodes)
  for (i in seq_along(run_patient_ids)) {
    m <- t(sapply(score_list, function(s) s[i, ]))                   # runs x nodes
    ranks <- t(apply(m, 1, function(v) rank(-v, ties.method = 'average')))
    for (r in seq_len(nrow(m))) rank_sum[[r]][i, ] <- ranks[r, ]
    rho_m <- suppressWarnings(cor(t(ranks)))
    top <- (t(apply(m, 1, function(v) rank(-v, ties.method = 'first'))) <= top_k) * 1
    inter <- top %*% t(top); jac <- inter / (2 * top_k - inter)
    patient_agreement[[i]] <- tibble(patient_id = run_patient_ids[i],
      real_real_rank_rho = mean(pair_mean(rho_m, is_init, is_init, TRUE), na.rm = TRUE),
      real_shuffled_rank_rho = mean(pair_mean(rho_m, is_init, is_perm, FALSE), na.rm = TRUE),
      real_real_topk_jaccard = mean(pair_mean(jac, is_init, is_init, TRUE)),
      real_shuffled_topk_jaccard = mean(pair_mean(jac, is_init, is_perm, FALSE)))
  }
  patient_agreement <- bind_rows(patient_agreement) %>%
    left_join(example_patients %>% select(patient_id, actual_exposure), by = 'patient_id')
  write_csv(patient_agreement, file.path(output_dir, 'Node_agreement_real_vs_shuffled_by_patient.csv'))

  # Cohort level: per run, median of each node's patient ranks (as in analyze.py).
  for (r in seq_len(nrow(agreement_runs))) cohort_rank_by_run[r, ] <- apply(rank_sum[[r]], 2, median)
  cohort_rho <- cor(t(cohort_rank_by_run))
  cohort_top <- (t(apply(cohort_rank_by_run, 1, function(v) rank(v, ties.method = 'first'))) <= top_k) * 1
  cohort_inter <- cohort_top %*% t(cohort_top); cohort_jac <- cohort_inter / (2 * top_k - cohort_inter)
  cohort_agreement <- bind_rows(
    tibble(comparison = 'Real vs real', measure = 'Rank Spearman', value = pair_mean(cohort_rho, is_init, is_init, TRUE)),
    tibble(comparison = 'Real vs shuffled', measure = 'Rank Spearman', value = pair_mean(cohort_rho, is_init, is_perm, FALSE)),
    tibble(comparison = 'Real vs real', measure = paste0('Top-', top_k, ' Jaccard'), value = pair_mean(cohort_jac, is_init, is_init, TRUE)),
    tibble(comparison = 'Real vs shuffled', measure = paste0('Top-', top_k, ' Jaccard'), value = pair_mean(cohort_jac, is_init, is_perm, FALSE)))
  write_csv(cohort_agreement, file.path(output_dir, 'Node_agreement_real_vs_shuffled_cohort.csv'))
  agreement_summary <- bind_rows(
    patient_agreement %>% summarise(level = 'patient (median over patients)',
      patients_real_real_higher_rho = sum(real_real_rank_rho > real_shuffled_rank_rho), n_patients = n(),
      real_real_rank_rho = median(real_real_rank_rho), real_shuffled_rank_rho = median(real_shuffled_rank_rho),
      real_real_topk_jaccard = median(real_real_topk_jaccard), real_shuffled_topk_jaccard = median(real_shuffled_topk_jaccard)),
    cohort_agreement %>% group_by(comparison, measure) %>% summarise(value = mean(value), .groups = 'drop') %>%
      mutate(key = paste(comparison, measure)) %>% select(key, value) %>%
      pivot_wider(names_from = key, values_from = value) %>%
      transmute(level = 'cohort (mean over run pairs)',
        real_real_rank_rho = `Real vs real Rank Spearman`, real_shuffled_rank_rho = `Real vs shuffled Rank Spearman`,
        real_real_topk_jaccard = .data[[paste0('Real vs real Top-', top_k, ' Jaccard')]],
        real_shuffled_topk_jaccard = .data[[paste0('Real vs shuffled Top-', top_k, ' Jaccard')]]))
  agreement_summary$top_k <- top_k
  write_csv(agreement_summary, file.path(output_dir, 'Node_agreement_summary.csv'))

  agreement_long <- patient_agreement %>%
    pivot_longer(-c(patient_id, actual_exposure), names_to = c('comparison', 'measure'),
                 names_pattern = '(real_real|real_shuffled)_(rank_rho|topk_jaccard)') %>%
    mutate(comparison = factor(recode(comparison, real_real = 'Real vs real', real_shuffled = 'Real vs shuffled'),
                               levels = c('Real vs real', 'Real vs shuffled')),
           measure = recode(measure, rank_rho = 'Node-rank Spearman', topk_jaccard = paste0('Top-', top_k, ' node Jaccard')))
  p_agree_patient <- ggplot(agreement_long, aes(comparison, value)) +
    geom_line(aes(group = patient_id), color = 'grey75', alpha = 0.6) +
    geom_boxplot(aes(fill = comparison), width = 0.4, outlier.shape = NA, alpha = 0.55) +
    geom_point(aes(color = comparison), size = 1.4, alpha = 0.8) +
    facet_wrap(~ measure, scales = 'free_y') +
    scale_fill_manual(values = c('#2886B5', '#A14792')) + scale_color_manual(values = c('#2886B5', '#A14792')) +
    labs(title = 'A. Patient-specific node rankings: real vs shuffled labels',
         subtitle = paste0(sum(is_init), ' initializations, ', sum(is_perm), ' permutations; lines join the same patient'),
         x = NULL, y = 'Mean agreement between run pairs') + theme(legend.position = 'none')
  p_agree_cohort <- ggplot(cohort_agreement, aes(comparison, value, fill = comparison)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(aes(color = comparison), width = 0.08, height = 0, size = 1, alpha = 0.6) +
    facet_wrap(~ measure, scales = 'free_y') +
    scale_fill_manual(values = c('#2886B5', '#A14792')) + scale_color_manual(values = c('#2886B5', '#A14792')) +
    labs(title = 'B. Cohort node rankings: real vs shuffled labels', x = NULL,
         y = 'Agreement per run pair',
         caption = 'Real vs real = seed-to-seed variability. Higher real-vs-real agreement means the\nhighlighted nodes depend on the exposure labels. Descriptive; no p-values.') +
    theme(legend.position = 'none')
  save_plot(p_agree_patient / p_agree_cohort, 'Node_agreement_real_vs_shuffled', 11, 10)
} else {
  message('Skipping real-vs-shuffled node agreement: needs >=2 initialization and >=1 permutation attributions.npz.')
}

# =============================================================================
# N3. Epoch robustness: compare with a second analysis run (e.g. 400 epochs)
# =============================================================================
if (nzchar(compare_analysis_dir)) {
  compare_config <- fromJSON(file.path(compare_analysis_dir, 'config.json'))
  this_name <- paste0(run_config$epochs, ' epochs'); other_name <- paste0(compare_config$epochs, ' epochs')
  if (this_name == other_name) { this_name <- 'this run'; other_name <- 'comparison run' }
  compare_nodes <- example_nodes %>% select(node_idx, node_label, rank_this = median_patient_rank) %>%
    inner_join(read_table(file.path(compare_analysis_dir, 'cohort_summary.csv')) %>%
                 select(node_idx, rank_other = median_patient_rank), by = 'node_idx')
  compare_edges <- example_edges %>% select(edge_id, rank_this = median_patient_rank) %>%
    inner_join(read_table(file.path(compare_analysis_dir, 'cohort_edge_scores.csv')) %>%
                 select(edge_id, rank_other = median_patient_rank), by = 'edge_id')
  top_overlap <- function(df, fraction) {
    k <- max(1, ceiling(fraction * nrow(df)))
    a <- df %>% arrange(rank_this) %>% slice_head(n = k) %>% pull(1)
    b <- df %>% arrange(rank_other) %>% slice_head(n = k) %>% pull(1)
    c(k = k, shared = length(intersect(a, b)), jaccard = length(intersect(a, b)) / length(union(a, b)))
  }
  gene_compare <- compare_nodes %>% filter(startsWith(node_label, 'GENE:'))
  robustness <- bind_rows(
    tibble(set = 'Gene nodes', spearman = rho(gene_compare$rank_this, gene_compare$rank_other),
           !!!as.list(top_overlap(gene_compare, node_fraction))),
    tibble(set = 'All nodes', spearman = rho(compare_nodes$rank_this, compare_nodes$rank_other),
           !!!as.list(top_overlap(compare_nodes, node_fraction))),
    tibble(set = 'Directed edges', spearman = rho(compare_edges$rank_this, compare_edges$rank_other),
           !!!as.list(top_overlap(compare_edges, edge_fraction)))) %>%
    mutate(this_run = this_name, comparison_run = other_name, top_fraction = c(node_fraction, node_fraction, edge_fraction))
  write_csv(robustness, file.path(output_dir, 'Epoch_robustness_summary.csv'))
  r_gene <- robustness %>% filter(set == 'Gene nodes')
  r_edge <- robustness %>% filter(set == 'Directed edges')
  p_rob_nodes <- ggplot(gene_compare, aes(rank_other, rank_this)) +
    geom_point(alpha = 0.25, size = 0.7, color = '#337A96') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
    labs(title = 'A. Gene-node cohort ranks', x = paste('Median patient rank,', other_name),
         y = paste('Median patient rank,', this_name),
         subtitle = paste0('Spearman = ', round(r_gene$spearman, 3), '; top ', round(node_fraction * 100),
                           '% shared = ', r_gene$shared, '/', r_gene$k))
  p_rob_edges <- ggplot(compare_edges, aes(rank_other, rank_this)) +
    geom_point(alpha = 0.08, size = 0.4, color = '#ED7658') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
    labs(title = 'B. Directed-edge cohort ranks', x = paste('Median patient rank,', other_name),
         y = paste('Median patient rank,', this_name),
         subtitle = paste0('Spearman = ', round(r_edge$spearman, 3), '; top ', round(edge_fraction * 100),
                           '% shared = ', r_edge$shared, '/', r_edge$k),
         caption = 'Lower rank = higher priority. Agreement indicates results are not sensitive to training length.')
  save_plot(p_rob_nodes + p_rob_edges, 'Epoch_robustness_node_edge_ranks', 12, 5.5)
}

# =============================================================================
# N4. Diagnostics from diagnose_exposure.py (optional DIAGNOSTICS_DIR)
# =============================================================================
if (nzchar(diagnostics_dir)) {
  diag_folders <- list.dirs(diagnostics_dir, recursive = FALSE)

  # Signal ceiling: cross-validated ridge vs shuffled-label null, per folder.
  ceiling_files <- file.path(diag_folders, 'ceiling_cv_ridge.csv')
  ceiling_files <- ceiling_files[file.exists(ceiling_files)]
  if (length(ceiling_files)) {
    ceil <- bind_rows(lapply(ceiling_files, function(f) read_table(f) %>% mutate(run = basename(dirname(f)))))
    write_csv(ceil, file.path(output_dir, 'Diagnostics_signal_ceiling.csv'))
    p_ceil <- ggplot(ceil, aes(y = feature_set)) +
      geom_errorbarh(aes(xmin = null_mean, xmax = null_95pct), height = 0.25, color = 'grey55') +
      geom_point(aes(x = null_mean), shape = 124, size = 4, color = 'grey45') +
      geom_point(aes(x = cv_within_fold_spearman), size = 3, color = '#2886B5') +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey60') +
      facet_wrap(~ run, ncol = 1) +
      labs(title = 'Exposure signal in the input data (no GNN)',
           x = 'Cross-validated Spearman (ridge)', y = NULL,
           caption = 'Blue = real labels; grey bar = shuffled-label mean to 95th percentile.\nNode-symmetric = what a model without node identity can see.')
    save_plot(p_ceil, 'Diagnostics_signal_ceiling', 9, 2.4 + 2 * length(ceiling_files))
  }

  # Known smoking markers, univariate.
  marker_files <- file.path(diag_folders, 'known_marker_associations.csv')
  marker_files <- marker_files[file.exists(marker_files)]
  if (length(marker_files)) {
    mk <- bind_rows(lapply(marker_files, function(f) read_table(f) %>% mutate(run = basename(dirname(f))))) %>%
      mutate(node_label = as.character(node_label))
    if (nrow(mk)) {
      p_mk <- ggplot(mk, aes(spearman, paste(node_label, modality), color = state)) +
        geom_vline(xintercept = 0, color = 'grey60') + geom_point(size = 2.5) +
        facet_wrap(~ run) + scale_color_manual(values = c(delta = '#ED7658', tumor = '#8361BB', normal = '#337A96')) +
        labs(title = 'Known smoking markers vs recorded exposure', x = 'Spearman with recorded exposure',
             y = NULL, color = 'Value', caption = 'Univariate, descriptive. Gene labels are model identifiers (e.g. Entrez IDs).')
      save_plot(p_mk, 'Diagnostics_known_markers', 10, 2 + 0.3 * n_distinct(paste(mk$node_label, mk$modality)))
    }
  }

  # Training trajectories: fitted rho vs epoch for real and shuffled labels.
  for (folder in diag_folders) {
    traj_file <- file.path(folder, 'trajectory.csv')
    if (!file.exists(traj_file)) next
    traj <- read_table(traj_file)
    traj_mean <- traj %>% group_by(config, labels, epoch) %>%
      summarise(rho = mean(fitted_spearman_vs_supplied), loss_ratio = mean(exposure_loss_vs_baseline), .groups = 'drop')
    gap <- traj_mean %>% select(config, labels, epoch, rho) %>% pivot_wider(names_from = labels, values_from = rho) %>%
      mutate(gap = observed - permuted)
    write_csv(gap, file.path(output_dir, paste0('Diagnostics_', basename(folder), '_gap_by_epoch.csv')))
    chosen_epoch <- run_config$epochs
    lab_colors <- c(observed = '#2886B5', permuted = '#A14792')
    p_traj <- ggplot(traj_mean, aes(epoch, rho, color = labels)) +
      geom_line(data = traj, aes(y = fitted_spearman_vs_supplied, group = interaction(labels, seed)),
                alpha = 0.25, linewidth = 0.3) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = chosen_epoch, linetype = 'dashed') +
      facet_wrap(~ config, nrow = 1) + scale_color_manual(values = lab_colors, labels = c('Real labels', 'Shuffled labels')) +
      labs(title = paste0('Training length: ', basename(folder)), y = 'In-sample Spearman', x = NULL, color = NULL) +
      theme(legend.position = 'top')
    p_gap <- ggplot(gap, aes(epoch, gap)) + geom_hline(yintercept = 0, color = 'grey60') +
      geom_line(color = 'black', linewidth = 0.8) + geom_vline(xintercept = chosen_epoch, linetype = 'dashed') +
      facet_wrap(~ config, nrow = 1) +
      labs(y = 'Real minus shuffled', x = 'Epoch',
           caption = paste0('Dashed line = epochs used in this analysis (', chosen_epoch,
                            '). Thin lines = seeds; thick = mean. In-sample diagnostics.'))
    p_loss <- ggplot(traj_mean, aes(epoch, loss_ratio, color = labels)) +
      geom_hline(yintercept = 1, linetype = 'dotted', color = 'grey50') + geom_line(linewidth = 0.8) +
      geom_vline(xintercept = chosen_epoch, linetype = 'dashed') +
      facet_wrap(~ config, nrow = 1) + scale_color_manual(values = lab_colors, labels = c('Real labels', 'Shuffled labels')) +
      labs(y = 'Exposure loss /\nconstant baseline', x = NULL, color = NULL) + theme(legend.position = 'none')
    n_cfg <- n_distinct(traj$config)
    save_plot(p_traj / p_loss / p_gap + plot_layout(heights = c(1.2, 0.8, 0.8)),
              paste0('Diagnostics_', basename(folder), '_training_length'), max(7, 3.2 * n_cfg), 9)
  }

  # Cross-seed node agreement from diagnose_exposure.py stages.
  agree_files <- file.path(diag_folders, 'node_agreement_across_seeds.csv')
  agree_files <- agree_files[file.exists(agree_files)]
  if (length(agree_files)) {
    ag <- bind_rows(lapply(agree_files, function(f) read_table(f) %>% mutate(run = basename(dirname(f)))))
    write_csv(ag, file.path(output_dir, 'Diagnostics_node_agreement_across_seeds.csv'))
  }
}

# =============================================================================
# Suppl. 6, 11, 12 and study limits
# =============================================================================
# Shared script exports training losses and initialization/subset/permutation
# confidence summaries. It cannot reconstruct missing per-epoch correlation logs.
# Patient score comparisons are descriptive. There is no independent validation set.
writeLines(c(
  'TCGA analysis uses the processed inputs in INPUT_DATA. No exposure values were reclassified;',
  'a filtered cohort (e.g. GNN_Input_Data_clean59) excludes patients rather than recoding them.',
  'Excluded patients: recorded zero exposure without a never-smoker status (see diagnostics patients_used.csv).',
  'Recorded zeros must not be equated with confirmed never-smoking without provenance.',
  'Original missingness, histology and other confounding remain unresolved when absent.',
  'An expectation of complete inputs does not establish completeness of this cohort.',
  'Available-cohort tables describe exported values; imputation can hide original missingness.',
  'Exposure enters only the loss. In-sample exposure fit is not independent evidence.',
  'Group comparisons and node/edge confidence measures are descriptive, without p-values.',
  'Edge scores describe fixed-network message suppression; they are not biological rewiring.',
  'Transition embeddings are exposure-trained and do not isolate exposure-independent cancer effects.',
  'Figure 3D pathway analysis is done separately; Supplement 9 is generated only when pathway_genes.csv is supplied.',
  'Real-vs-shuffled node agreement and epoch-robustness plots are descriptive comparisons, not tests.',
  'Score emphasises tumor-minus-normal changes; exposure effects shared by both tissues are out of scope.',
  'No adjusted/confounder-independent or out-of-sample predictive claims are supported by these plots.',
  'Supplement 10 uses candidate CpG-gene channels; regulatory validity is not established.'
),file.path(output_dir,'TCGA_limitations.txt'))
message('TCGA-specific plots completed: ', normalizePath(output_dir))
