# =============================================================================
# Held-out evaluation and ranking stability across splits and seeds
# =============================================================================
# Reads ExposureGAT model outputs only (no refitting; model scores unchanged).
#
# Subset fits train on a random 80% of patients; the remaining 20% are scored by
# the fitted model without having been seen (scalers and cohort reference are
# also fitted on the training patients only). This script reports:
#   1. Exposure-rank agreement per split: training patients vs held-out patients.
#   2. Pooled held-out agreement: each patient's mean score over the fits that
#      excluded them, against recorded exposure (with a label-shuffle reference).
#   3. Fit agreement across random seeds (initialization fits).
#   4. Patient-level ranking stability: how often each patient's consensus top
#      nodes stay in their top list across seeds, and across splits when that
#      patient was included vs held out.
#   5. Cohort ranking stability: pairwise agreement of cohort node rankings
#      across seeds and across splits.
# Writes tables and one figure.
#
#   Rscript held_out_summary.R ANALYSIS_DIR OUTPUT_DIR
# Example:
#   Rscript held_out_summary.R "$V/TCGA_clean59_e600/analysis" "$V/TCGA_clean59_e600/held_out"
# install.packages(c('dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite'))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript held_out_summary.R ANALYSIS_DIR OUTPUT_DIR')
analysis_dir <- args[1]; output_dir <- args[2]
if (!dir.exists(analysis_dir)) stop('Folder not found: ', analysis_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) + theme(plot.title = element_text(face = 'bold')))
read_table <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)
set.seed(42) # label-shuffle reference only
config <- fromJSON(file.path(analysis_dir, 'config.json'))

rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = 'spearman')
}
fmt <- function(x, d = 2) formatC(x, format = 'f', digits = d)
med_iqr <- function(x) paste0(fmt(median(x, na.rm = TRUE)), ' (IQR ', fmt(quantile(x, .25, na.rm = TRUE)),
                              ' to ', fmt(quantile(x, .75, na.rm = TRUE)), ')')

read_npy_from_npz <- function(npz, key) {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(npz, files = paste0(key, '.npy'), exdir = tmp)
  path <- file.path(tmp, paste0(key, '.npy'))
  bytes <- readBin(path, 'raw', file.info(path)$size)
  if (!identical(bytes[1:6], as.raw(c(0x93, 0x4E, 0x55, 0x4D, 0x50, 0x59)))) stop('Not a .npy file: ', key)
  if (as.integer(bytes[7]) == 1) { hlen <- sum(as.integer(bytes[9:10]) * 256^(0:1)); start <- 11 }
  else { hlen <- sum(as.integer(bytes[9:12]) * 256^(0:3)); start <- 13 }
  header <- rawToChar(bytes[start:(start + hlen - 1)])
  descr <- sub(".*'descr': *'([^']+)'.*", '\\1', header)
  shape <- as.integer(strsplit(gsub('[^0-9,]', '', sub(".*'shape': *\\(([^)]*)\\).*", '\\1', header)), ',')[[1]])
  shape <- shape[!is.na(shape)]
  size <- switch(descr, '<f4' = 4, '<f8' = 8, stop('Unsupported dtype ', descr))
  values <- readBin(bytes[(start + hlen):length(bytes)], 'numeric', n = prod(shape), size = size, endian = 'little')
  if (length(shape) == 2) matrix(values, nrow = shape[1], ncol = shape[2], byrow = TRUE) else values
}

# =============================================================================
# 1-3. Exposure-rank agreement: per split, pooled held-out, per seed
# =============================================================================
runs_scores <- read_table(file.path(analysis_dir, 'patient_scores_by_run.csv')) %>%
  mutate(patient_id = as.character(patient_id), included_in_fit = as.logical(included_in_fit))
subset <- runs_scores %>% filter(kind == 'subset')
init <- runs_scores %>% filter(kind == 'initialization')
if (!nrow(subset)) stop('No subset fits found in ', analysis_dir)

per_split <- subset %>% group_by(`repeat`) %>% summarise(
  n_training = sum(included_in_fit), n_held_out = sum(!included_in_fit),
  training_rho = rho(predicted_exposure_rank_score[included_in_fit], supplied_exposure[included_in_fit]),
  held_out_rho = rho(predicted_exposure_rank_score[!included_in_fit], supplied_exposure[!included_in_fit]),
  .groups = 'drop')
write_csv(per_split, file.path(output_dir, 'Held_out_by_split.csv'))

pooled <- subset %>% filter(!included_in_fit) %>% group_by(patient_id) %>%
  summarise(mean_held_out_score = mean(predicted_exposure_rank_score), times_held_out = n(),
            exposure = first(supplied_exposure), .groups = 'drop')
pooled_rho <- rho(pooled$mean_held_out_score, pooled$exposure)
shuffle_ref <- replicate(2000, rho(pooled$mean_held_out_score, sample(pooled$exposure)))
pooled_shuffle_fraction <- mean(shuffle_ref >= pooled_rho, na.rm = TRUE)
write_csv(pooled, file.path(output_dir, 'Held_out_pooled_by_patient.csv'))

per_seed <- init %>% group_by(`repeat`) %>%
  summarise(fitted_rho = rho(predicted_exposure_rank_score, supplied_exposure), .groups = 'drop')
write_csv(per_seed, file.path(output_dir, 'Fit_by_seed.csv'))

# =============================================================================
# 4. Patient-level ranking stability of each patient's consensus top nodes
# =============================================================================
k <- config$top_k
patient_nodes <- read_table(file.path(analysis_dir, 'patient_node_scores.csv')) %>%
  mutate(patient_id = as.character(patient_id))
patient_stability <- patient_nodes %>% group_by(patient_id) %>%
  arrange(median_rank, node_idx, .by_group = TRUE) %>% slice_head(n = k) %>%
  summarise(across_seeds = mean(initialization_topk_frequency),
            across_splits_included = mean(subset_included_frequency, na.rm = TRUE),
            across_splits_held_out = mean(subset_excluded_frequency, na.rm = TRUE), .groups = 'drop')
write_csv(patient_stability, file.path(output_dir, 'Patient_top_node_stability.csv'))
# Reference: chance that a node is in a random top-k list.
n_nodes <- n_distinct(patient_nodes$node_idx); chance_topk <- k / n_nodes

# =============================================================================
# 5. Cohort ranking stability across seeds and across splits
# =============================================================================
run_table <- read_table(file.path(analysis_dir, 'runs.csv')) %>% filter(kind %in% c('initialization', 'subset')) %>%
  mutate(folder = file.path(analysis_dir, paste0(kind, '_', .data[['repeat']])))
cohort_rank <- t(sapply(run_table$folder, function(f) {
  s <- read_npy_from_npz(file.path(f, 'attributions.npz'), 'score')
  r <- t(apply(s, 1, function(v) rank(-v, ties.method = 'average')))   # patient x node ranks
  apply(r, 2, median)
}))
pair_stats <- function(idx) {
  if (length(idx) < 2) return(tibble())
  combos <- combn(idx, 2)
  bind_rows(lapply(seq_len(ncol(combos)), function(j) {
    a <- cohort_rank[combos[1, j], ]; b <- cohort_rank[combos[2, j], ]
    ta <- order(a)[seq_len(k)]; tb <- order(b)[seq_len(k)]
    tibble(spearman = cor(a, b, method = 'spearman'),
           topk_jaccard = length(intersect(ta, tb)) / length(union(ta, tb)))
  }))
}
cohort_stability <- bind_rows(
  pair_stats(which(run_table$kind == 'initialization')) %>% mutate(comparison = 'Across seeds'),
  pair_stats(which(run_table$kind == 'subset')) %>% mutate(comparison = 'Across 80% splits'))
write_csv(cohort_stability, file.path(output_dir, 'Cohort_ranking_stability.csv'))

# =============================================================================
# Summary table and figure
# =============================================================================
summary <- tibble(
  measure = c('Fitted exposure-rank Spearman, across seeds',
              'Training-patient Spearman, across 80% splits',
              'Held-out-patient Spearman, across 80% splits',
              'Pooled held-out Spearman (all patients ever held out)',
              'Share of label shuffles >= pooled held-out Spearman',
              'Patient top-node retention across seeds',
              'Patient top-node retention across splits, patient included',
              'Patient top-node retention across splits, patient held out',
              'Chance retention for a random top list',
              'Cohort ranking Spearman across seeds',
              'Cohort ranking Spearman across splits',
              'Cohort top-k Jaccard across seeds',
              'Cohort top-k Jaccard across splits'),
  value = c(med_iqr(per_seed$fitted_rho), med_iqr(per_split$training_rho), med_iqr(per_split$held_out_rho),
            paste0(fmt(pooled_rho), ' (n = ', nrow(pooled), ' patients)'), fmt(pooled_shuffle_fraction, 3),
            med_iqr(patient_stability$across_seeds), med_iqr(patient_stability$across_splits_included),
            med_iqr(patient_stability$across_splits_held_out), fmt(chance_topk, 3),
            med_iqr(cohort_stability$spearman[cohort_stability$comparison == 'Across seeds']),
            med_iqr(cohort_stability$spearman[cohort_stability$comparison == 'Across 80% splits']),
            med_iqr(cohort_stability$topk_jaccard[cohort_stability$comparison == 'Across seeds']),
            med_iqr(cohort_stability$topk_jaccard[cohort_stability$comparison == 'Across 80% splits'])))
write_csv(summary, file.path(output_dir, 'Held_out_summary.csv'))
print(as.data.frame(summary), right = FALSE)

fit_long <- bind_rows(per_seed %>% transmute(group = 'Training\n(all patients, seeds)', value = fitted_rho),
                      per_split %>% transmute(group = 'Training\n(80% splits)', value = training_rho),
                      per_split %>% transmute(group = 'Held out\n(20% splits)', value = held_out_rho)) %>%
  mutate(group = factor(group, levels = c('Training\n(all patients, seeds)', 'Training\n(80% splits)', 'Held out\n(20% splits)')))
p1 <- ggplot(fit_long, aes(group, value)) + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_boxplot(outlier.shape = NA, width = 0.45, fill = '#DCE9F2') +
  geom_jitter(width = 0.08, height = 0, size = 1.8, alpha = 0.8, color = '#2886B5') +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(title = 'A. Exposure-rank agreement', x = NULL, y = 'Spearman (score vs recorded exposure)',
       subtitle = paste0('Pooled held-out Spearman = ', fmt(pooled_rho), ' (n = ', nrow(pooled), ')'))
stab_long <- patient_stability %>% pivot_longer(-patient_id, names_to = 'group', values_to = 'value') %>%
  filter(is.finite(value)) %>%
  mutate(group = factor(recode(group, across_seeds = 'Across seeds', across_splits_included = 'Splits,\nincluded',
                               across_splits_held_out = 'Splits,\nheld out'),
                        levels = c('Across seeds', 'Splits,\nincluded', 'Splits,\nheld out')))
p2 <- ggplot(stab_long, aes(group, value)) + geom_hline(yintercept = chance_topk, linetype = 'dashed', color = 'grey50') +
  geom_boxplot(outlier.shape = NA, width = 0.45, fill = '#F2E3DC') +
  geom_jitter(width = 0.08, height = 0, size = 1.2, alpha = 0.6, color = '#C17A26') +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = paste0('B. Patient top-', k, ' node retention'), x = NULL,
       y = 'Fraction of fits keeping the node in top list', subtitle = 'Dashed = chance')
p3 <- ggplot(cohort_stability, aes(comparison, spearman)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, fill = '#E8E0F0') +
  geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.5, color = '#8361BB') +
  coord_cartesian(ylim = c(-0.2, 1)) +
  labs(title = 'C. Cohort ranking agreement', x = NULL, y = 'Spearman between run pairs')
fig <- p1 + p2 + p3 + plot_annotation(caption = paste0(config$epochs, ' epochs; ', nrow(per_seed), ' seeds; ',
  nrow(per_split), ' splits of 80%/20%. Held-out patients were not used to fit model parameters, scalers or the cohort reference.'))
ggsave(file.path(output_dir, 'Held_out_and_stability.pdf'), fig, width = 15, height = 5.5)
ggsave(file.path(output_dir, 'Held_out_and_stability.png'), fig, width = 15, height = 5.5, dpi = 200)


message('Saved to ', normalizePath(output_dir))
