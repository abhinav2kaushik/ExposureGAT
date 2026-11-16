# =============================================================================
# Cohort analysis (TCGA): Figure 3A and Supplementary Figures 6-8
# =============================================================================
#   Rscript cohort_figures.R DATA_DIR MODEL_DIR FIGURE_DIR
# Figure 3A: in-sample fit of the patient score to the standardized exposure rank.
# Supplementary Figure 6: training losses. Supplementary Figure 7: distributions of gene-node
# and edge priorities. Supplementary Figure 8: node priority vs incoming edge priority per patient.
# Figures 3B-C come from src/scoring/weighted_node_score.py.
# Reads model outputs only; no fitting and no change to model scores.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('Usage: Rscript cohort_figures.R DATA_DIR MODEL_DIR FIGURE_DIR')
data_dir <- args[1]; model_dir <- args[2]; fig_dir <- args[3]
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, '.pdf')), plot, width = width, height = height)
  ggsave(file.path(fig_dir, paste0(name, '.png')), plot, width = width, height = height, dpi = 300)
}
rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = 'spearman')
}
theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
                                           plot.title = element_text(face = 'bold')))
kind_colors <- c(initialization = '#2886B5', subset = '#C17A26', permutation = '#A14792')
kind_labels <- c(initialization = 'Initialization', subset = 'Subset', permutation = 'Permutation')
config <- fromJSON(file.path(model_dir, 'config.json'))

# =============================================================================
# Figure 3A: fitted vs recorded standardized exposure rank (median across initializations)
# =============================================================================
meta <- read_table(file.path(data_dir, 'metadata_samples.csv')) %>%
  transmute(patient_id = as.character(patient_id), exposure = .data[[config$exposure_col]]) %>% distinct()
patients <- read_table(file.path(model_dir, 'patient_scores.csv')) %>% mutate(patient_id = as.character(patient_id)) %>%
  left_join(meta, by = 'patient_id') %>%
  mutate(q = (rank(exposure, ties.method = 'average') - 0.5) / n(),
         recorded_rank = (q - mean(q)) / sqrt(mean((q - mean(q))^2)),   # population SD, as in training
         fitted_rank = predicted_exposure_rank_score_median)
write_csv(select(patients, patient_id, exposure, recorded_rank, fitted_rank), file.path(fig_dir, 'Figure_3A_data.csv'))
r3a <- rho(patients$fitted_rank, patients$recorded_rank)
lim <- range(c(patients$recorded_rank, patients$fitted_rank))
p3a <- ggplot(patients, aes(recorded_rank, fitted_rank)) +
  geom_abline(linetype = 'dashed', color = 'grey50') + geom_point(alpha = 0.7, color = '#8361BB') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black') +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = 'Fitted vs recorded exposure rank', subtitle = paste0('Spearman rho = ', sprintf('%.3f', r3a)),
       x = 'Recorded standardized exposure rank', y = 'Fitted standardized exposure-rank score')
save_plot(p3a, 'Figure_3A', 5, 5)

# =============================================================================
# Supplementary Figure 6: training losses of all fits
# =============================================================================
runs <- read_table(file.path(model_dir, 'runs.csv'))
h <- bind_rows(lapply(seq_len(nrow(runs)), function(j)
  read_table(file.path(model_dir, paste0(runs$kind[j], '_', runs[['repeat']][j]), 'history.csv')) %>%
    mutate(kind = runs$kind[j], repetition = runs[['repeat']][j]))) %>%
  pivot_longer(c(exposure_loss, transition_loss), names_to = 'objective', values_to = 'value') %>%
  mutate(objective = factor(objective, levels = c('exposure_loss', 'transition_loss'),
                            labels = c('Exposure-rank loss', 'Transition-reconstruction loss')))
h_mean <- h %>% group_by(kind, objective, epoch) %>% summarise(value = mean(value), .groups = 'drop')
s6 <- ggplot(h, aes(epoch, value, color = kind)) +
  geom_line(aes(group = interaction(kind, repetition)), alpha = 0.15, linewidth = 0.3) +
  geom_line(data = h_mean, linewidth = 0.9) + facet_wrap(~objective, ncol = 1, scales = 'free_y') +
  scale_color_manual(values = kind_colors, labels = kind_labels) +
  labs(title = 'Model optimization', x = 'Epoch', y = 'Training loss', color = NULL)
save_plot(s6, 'Supplementary_Figure_06', 8, 8)

# =============================================================================
# Supplementary Figure 7: distributions of cohort gene-node and edge priorities
# =============================================================================
# Normalized priority = (N - median patient rank) / (N - 1); the dashed line marks the top-10% cutoff.
nodes <- read_table(file.path(model_dir, 'cohort_summary.csv'))
edges <- read_table(file.path(model_dir, 'cohort_edge_scores.csv'))
normalize <- function(r) (length(r) - r) / max(length(r) - 1, 1)
genes <- nodes %>% mutate(priority = normalize(median_patient_rank)) %>% filter(startsWith(node_label, 'GENE:'))  # normalized over all nodes
edges <- edges %>% mutate(priority = normalize(median_patient_rank))
distribution <- function(values, title, x_label) {
  cutoff <- sort(values, decreasing = TRUE)[ceiling(0.10 * length(values))]
  ggplot(tibble(x = values), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = 'grey75', color = 'grey30') +
    geom_density(linewidth = 0.7) + geom_vline(xintercept = cutoff, linetype = 'dashed', color = 'red') +
    annotate('text', x = Inf, y = Inf, label = paste0('Top 10% cutoff = ', signif(cutoff, 3)),
             hjust = 1.05, vjust = 1.5, size = 3.5, color = 'red') +
    labs(title = title, x = x_label, y = 'Density')
}
s7 <- distribution(genes$priority, 'Gene-node priority', 'Normalized node priority') /
  distribution(edges$priority, 'Edge priority', 'Normalized edge priority') + plot_annotation(tag_levels = 'A')
save_plot(s7, 'Supplementary_Figure_07', 8, 9)

# =============================================================================
# Supplementary Figure 8: node priority vs incoming edge priority, one point per patient
# =============================================================================
# Per patient: take the top 10% edges by edge priority, sum their priorities at each target node,
# take the top 10% nodes by this incoming strength, and compare their mean incoming strength with
# their median node priority. Both summaries come from the same fitted model.
patient_nodes <- read_table(file.path(model_dir, 'patient_node_scores.csv')) %>% mutate(patient_id = as.character(patient_id))
patient_edges <- read_csv(file.path(model_dir, 'patient_edge_scores.csv.gz'), show_col_types = FALSE, progress = FALSE,
                          col_types = cols_only(patient_id = col_character(), edge_id = col_integer(),
                                                median_rank = col_double(), median_priority = col_double()))
strength <- patient_edges %>% group_by(patient_id) %>% arrange(median_rank, edge_id, .by_group = TRUE) %>%
  slice_head(n = ceiling(0.10 * nrow(edges))) %>% ungroup() %>%
  left_join(select(edges, edge_id, target_idx), by = 'edge_id') %>%
  group_by(patient_id, node_idx = target_idx) %>% summarise(strength = sum(median_priority), .groups = 'drop')
s8_data <- patient_nodes %>% select(patient_id, node_idx, median_priority) %>%
  left_join(strength, by = c('patient_id', 'node_idx')) %>% mutate(strength = replace_na(strength, 0)) %>%
  group_by(patient_id) %>% arrange(desc(strength), node_idx, .by_group = TRUE) %>%
  slice_head(n = ceiling(0.10 * nrow(nodes))) %>%
  summarise(edge_priority = mean(strength), node_priority = median(median_priority), .groups = 'drop')
write_csv(s8_data, file.path(fig_dir, 'Supplementary_Figure_08_data.csv'))
s8 <- ggplot(s8_data, aes(edge_priority, node_priority)) +
  geom_point(size = 2, alpha = 0.7, color = 'grey40') +
  geom_smooth(method = 'lm', formula = y ~ x, color = '#4278FF', fill = 'grey70') +
  labs(title = 'Edge vs node priority',
       subtitle = paste0('Spearman rho = ', sprintf('%.3f', rho(s8_data$edge_priority, s8_data$node_priority))),
       x = 'Edge priority at top nodes (mean incoming, top 10% edges)', y = 'Node priority (median, top 10% nodes)')
save_plot(s8, 'Supplementary_Figure_08', 7, 6)

message('Cohort figures written to ', normalizePath(fig_dir))
