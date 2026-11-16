# =============================================================================
# Simulation benchmark: Figure 2B-D and Supplementary Figures 1-5
# =============================================================================
#   Rscript simulation_figures.R MANIFEST.csv FIGURE_DIR
# MANIFEST.csv has one row per independent simulated dataset (columns seed, data_dir,
# analysis_dir). The first row is the illustrated dataset (Figure 2B, Supplementary
# Figures 2-4); it is listed first by design, never chosen by performance.
# Figure 2A is the design schematic and is not generated here.
# Reads model outputs only; no fitting and no change to model scores.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(jsonlite); library(ggpubr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript simulation_figures.R MANIFEST.csv FIGURE_DIR')
manifest <- read_csv(args[1], show_col_types = FALSE) %>% mutate(seed = as.character(seed))
fig_dir <- args[2]
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(c('seed', 'data_dir', 'analysis_dir') %in% names(manifest)), !anyDuplicated(manifest$seed))
example_seed <- manifest$seed[1]
set.seed(42) # plotting jitter only

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
set_colors <- c(A = '#E8734E', C = '#704B95', B = '#8A8A8A')
set_labels <- c(A = 'Set-A', C = 'Set-C', B = 'Set-B')
kind_colors <- c(initialization = '#2886B5', subset = '#C17A26', permutation = '#A14792')
kind_labels <- c(initialization = 'Initialization', subset = 'Subset', permutation = 'Permutation')

# =============================================================================
# Read every simulated dataset
# =============================================================================
patients <- list(); nodes <- list(); histories <- list(); degrees <- list(); audit <- list()
for (i in seq_len(nrow(manifest))) {
  d <- manifest$data_dir[i]; a <- manifest$analysis_dir[i]; id <- manifest$seed[i]
  config <- fromJSON(file.path(a, 'config.json'))
  meta <- read_table(file.path(d, 'metadata_samples.csv')) %>%
    transmute(patient_id = as.character(patient_id), exposure = .data[[config$exposure_col]]) %>% distinct()
  truth <- read_table(file.path(d, 'sim_ground_truth_nodes.csv')) %>% select(node_idx, set_membership)

  p <- read_table(file.path(a, 'patient_scores.csv')) %>% mutate(patient_id = as.character(patient_id)) %>%
    left_join(meta, by = 'patient_id') %>%
    mutate(seed = id, normalized_exposure = (rank(exposure, ties.method = 'average') - 0.5) / n(),
           transition_magnitude = transition_score_median)
  patients[[i]] <- p

  nodes[[i]] <- read_table(file.path(a, 'cohort_summary.csv')) %>% left_join(truth, by = 'node_idx') %>%
    left_join(read_table(file.path(a, 'cohort_node_transition_scores.csv')), by = c('node_idx', 'node_label')) %>%
    mutate(seed = id, node_priority = (n() - median_patient_rank) / max(n() - 1, 1))

  runs <- read_table(file.path(a, 'runs.csv'))
  for (j in seq_len(nrow(runs))) {
    histories[[length(histories) + 1]] <- read_table(file.path(a, paste0(runs$kind[j], '_', runs[['repeat']][j]), 'history.csv')) %>%
      mutate(seed = id, kind = runs$kind[j], repetition = runs[['repeat']][j])
  }

  # Unique undirected degree of every node in the graph supplied to the model.
  und <- read_table(file.path(d, 'edges_ppi.csv')) %>% transmute(u = pmin(from_idx, to_idx), v = pmax(from_idx, to_idx)) %>%
    filter(u != v) %>% distinct()
  deg <- bind_rows(transmute(und, node_idx = u), transmute(und, node_idx = v)) %>% count(node_idx, name = 'degree')
  degrees[[i]] <- truth %>% left_join(deg, by = 'node_idx') %>% mutate(degree = replace_na(degree, 0L), seed = id)

  # Total observed paired change per modality, and after the loader transformation.
  samples <- read_table(file.path(d, 'metadata_samples.csv')) %>% mutate(across(c(sample_id, patient_id), as.character))
  pairs <- inner_join(samples %>% filter(tissue_type == 'NT') %>% select(patient_id, normal = sample_id),
                      samples %>% filter(tissue_type == 'TP') %>% select(patient_id, tumor = sample_id), by = 'patient_id')
  for (m in c('rna', 'meth', 'cnv')) {
    x <- read_table(file.path(d, paste0('X_', m, '_gene_features.csv'))); names(x)[1] <- 'sample_id'
    xm <- as.matrix(x[, -1]); delta <- xm[match(pairs$tumor, x$sample_id), ] - xm[match(pairs$normal, x$sample_id), ]
    audit[[length(audit) + 1]] <- tibble(seed = id, patient_id = pairs$patient_id, measure = m,
                                         magnitude = rowSums(delta^2)) %>% left_join(meta, by = 'patient_id')
  }
  audit[[length(audit) + 1]] <- p %>% transmute(seed, patient_id, measure = 'loader',
                                                magnitude = scaled_input_transition_rms_median, exposure)
}
patients <- bind_rows(patients); nodes <- bind_rows(nodes); histories <- bind_rows(histories)
degrees <- bind_rows(degrees); audit <- bind_rows(audit)
example_patients <- filter(patients, seed == example_seed)
example_nodes <- filter(nodes, seed == example_seed) %>% mutate(set = factor(set_membership, levels = c('A', 'C', 'B')))

# =============================================================================
# Figure 2B: patient-level paired-transition magnitude vs exposure (illustrated dataset)
# =============================================================================
r2b <- rho(example_patients$transition_magnitude, example_patients$normalized_exposure)
p2b <- ggplot(example_patients, aes(normalized_exposure, transition_magnitude)) +
  geom_point(alpha = 0.7, color = '#337A96') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black') +
  labs(title = 'Paired-transition magnitude vs exposure', subtitle = paste0('Spearman rho = ', round(r2b, 2)),
       x = 'Normalized exposure rank', y = 'Mean paired embedding change')
save_plot(p2b, 'Figure_2B', 5, 4)

# =============================================================================
# Figure 2C: node-level ROC curves, averaged across datasets
# =============================================================================
comparisons <- list('A vs C' = c('A', 'C'), 'A vs B' = c('A', 'B'), 'C vs B' = c('C', 'B'))
grid <- seq(0, 1, length.out = 501); roc <- list(); aucs <- list()
for (id in manifest$seed) for (cmp in names(comparisons)) {
  cl <- comparisons[[cmp]]
  x <- nodes %>% filter(seed == id, set_membership %in% cl) %>% mutate(pos = set_membership == cl[1])
  g <- x %>% group_by(node_priority) %>% summarise(tp = sum(pos), fp = sum(!pos), .groups = 'drop') %>% arrange(desc(node_priority))
  fpr <- c(0, cumsum(g$fp) / sum(!x$pos)); tpr <- c(0, cumsum(g$tp) / sum(x$pos))
  aucs[[length(aucs) + 1]] <- tibble(seed = id, comparison = cmp, AUROC = sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
  k <- tibble(fpr, tpr) %>% group_by(fpr) %>% summarise(tpr = max(tpr), .groups = 'drop')
  roc[[length(roc) + 1]] <- tibble(seed = id, comparison = cmp, fpr = grid, tpr = approx(k$fpr, k$tpr, xout = grid, rule = 2)$y)
}
auc_summary <- bind_rows(aucs) %>% group_by(comparison) %>% summarise(mean = mean(AUROC), sd = sd(AUROC), .groups = 'drop')
roc_mean <- bind_rows(roc) %>% group_by(comparison, fpr) %>% summarise(tpr = mean(tpr), .groups = 'drop')
roc_mean <- bind_rows(tibble(comparison = names(comparisons), fpr = 0, tpr = 0), roc_mean)
write_csv(bind_rows(aucs), file.path(fig_dir, 'Figure_2C_AUROC_by_dataset.csv'))
labels_2c <- setNames(paste0(auc_summary$comparison, '  AUROC: ', sprintf('%.2f', auc_summary$mean)), auc_summary$comparison)
p2c <- ggplot(roc_mean, aes(fpr, tpr, color = comparison)) +
  geom_abline(linetype = 'dashed') + geom_line(linewidth = 1.2) +
  scale_color_manual(values = c('A vs C' = '#227C9D', 'A vs B' = '#E8734E', 'C vs B' = '#704B95'),
                     breaks = names(comparisons), labels = labels_2c[names(comparisons)]) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = 'Node-level ROC', subtitle = 'Curves averaged across datasets',
       x = 'False positive rate', y = 'True positive rate', color = NULL) +
  theme(legend.position = c(.64, .15), legend.background = element_rect(fill = 'white'))
save_plot(p2c, 'Figure_2C', 6.5, 6)

# =============================================================================
# Figure 2D: composition of the top-50 nodes vs random-selection expectation
# =============================================================================
top_k <- 50
counts <- nodes %>% group_by(seed) %>% arrange(median_patient_rank, node_idx, .by_group = TRUE) %>%
  mutate(selected = row_number() <= top_k) %>%
  group_by(seed, set_membership) %>% summarise(selected = sum(selected), set_size = n(), .groups = 'drop') %>%
  group_by(seed) %>% mutate(expected = sum(selected) * set_size / sum(set_size)) %>% ungroup() %>%
  mutate(set = factor(set_membership, levels = c('A', 'C', 'B')))
write_csv(counts, file.path(fig_dir, 'Figure_2D_top50_counts.csv'))
means <- counts %>% group_by(set) %>% summarise(selected = mean(selected), expected = mean(expected), .groups = 'drop')
p2d <- ggplot(means, aes(set, selected, fill = set)) + geom_col(alpha = 0.65) +
  geom_point(aes(y = expected), shape = 95, size = 14) +
  geom_point(data = counts, size = 2, position = position_jitter(width = 0.25, height = 0)) +
  scale_fill_manual(values = setNames(set_colors, c('A', 'C', 'B'))) + scale_x_discrete(labels = set_labels) +
  labs(title = paste0('Top-', top_k, ' nodes'), x = NULL, y = 'Number selected') + theme(legend.position = 'none')
save_plot(p2d, 'Figure_2D', 3.6, 4)

# =============================================================================
# Supplementary Figure 1: total observed change vs exposure (A) and input connectivity (B)
# =============================================================================
measure_labels <- c(rna = 'RNA (L2 squared)', meth = 'DNA methylation (L2 squared)',
                    cnv = 'CNV (L2 squared)', loader = 'Loader-transformed RMS')
s1a_data <- audit %>% group_by(seed, measure) %>% summarise(rho = rho(magnitude, exposure), .groups = 'drop') %>%
  mutate(measure = factor(measure_labels[measure], levels = measure_labels))
write_csv(s1a_data, file.path(fig_dir, 'Supplementary_Figure_01A_data.csv'))
s1a <- ggplot(s1a_data, aes(measure, rho)) + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_jitter(width = 0.08, height = 0, size = 2) + coord_cartesian(ylim = c(-1, 1)) +
  labs(title = 'Total observed change vs exposure', x = NULL, y = 'Spearman correlation per dataset') +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
s1b <- ggplot(degrees, aes(factor(set_membership, levels = c('A', 'C', 'B'), labels = set_labels), degree)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 1) +
  labs(title = 'Actual input connectivity', x = NULL, y = 'Unique undirected degree')
save_plot(s1a / s1b + plot_annotation(tag_levels = 'A'), 'Supplementary_Figure_01', 7, 9)

# =============================================================================
# Supplementary Figure 2: training losses (illustrated dataset)
# =============================================================================
h <- histories %>% filter(seed == example_seed) %>%
  pivot_longer(c(exposure_loss, transition_loss), names_to = 'objective', values_to = 'value') %>%
  mutate(objective = factor(objective, levels = c('exposure_loss', 'transition_loss'),
                            labels = c('Exposure-rank loss', 'Transition-reconstruction loss')))
h_mean <- h %>% group_by(kind, objective, epoch) %>% summarise(value = mean(value), .groups = 'drop')
s2 <- ggplot(h, aes(epoch, value, color = kind)) +
  geom_line(aes(group = interaction(kind, repetition)), alpha = 0.15, linewidth = 0.3) +
  geom_line(data = h_mean, linewidth = 0.9) + facet_wrap(~objective, ncol = 1, scales = 'free_y') +
  scale_color_manual(values = kind_colors, labels = kind_labels) +
  labs(title = 'Model optimization', x = 'Epoch', y = 'Training loss', color = NULL)
save_plot(s2, 'Supplementary_Figure_02', 7, 7)

# =============================================================================
# Supplementary Figures 3-4: node scores by simulated node class (illustrated dataset)
# =============================================================================
class_violin <- function(column, title, y_label) {
  pairs <- list(c('A', 'C'), c('A', 'B'), c('C', 'B'))
  y_max <- max(example_nodes[[column]], na.rm = TRUE)
  ggplot(example_nodes, aes(set, .data[[column]], fill = set)) +
    geom_violin(alpha = 0.1, trim = FALSE) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.05, size = 1.2, alpha = 0.6) +
    scale_fill_manual(values = setNames(set_colors, c('A', 'C', 'B'))) + scale_x_discrete(labels = set_labels) +
    stat_compare_means(comparisons = pairs, method = 'wilcox.test', label = 'p.format', tip.length = 0.01) +
    stat_compare_means(method = 'kruskal.test', label.x = 0.6, label.y = y_max * 1.3) +
    labs(title = title, x = NULL, y = y_label) + theme(legend.position = 'none')
}
save_plot(class_violin('node_priority', 'Cohort node-priority score by node class', 'Cohort node-priority score'),
          'Supplementary_Figure_03', 6.5, 7)
save_plot(class_violin('median_patient_transition_score', 'Paired-transition magnitude by node class', 'Paired-transition magnitude'),
          'Supplementary_Figure_04', 6.5, 7)

# =============================================================================
# Supplementary Figure 5: exposure association of edge-level attention changes
# =============================================================================
# Edge score = |Spearman| between the patient's signed gated-attention change of an
# undirected connection (two directions averaged) and exposure. Edge classes refer to the
# classes of the two endpoint nodes; the simulation plants no edge-specific changes.
edge_scores <- list(); degree_rows <- list()
for (i in seq_len(nrow(manifest))) {
  d <- manifest$data_dir[i]; a <- manifest$analysis_dir[i]; id <- manifest$seed[i]
  config <- fromJSON(file.path(a, 'config.json'))
  truth <- read_table(file.path(d, 'sim_ground_truth_nodes.csv')) %>% select(node_idx, set_membership)
  meta <- read_table(file.path(d, 'metadata_samples.csv')) %>%
    transmute(patient_id = as.character(patient_id), exposure = .data[[config$exposure_col]]) %>% distinct()
  e <- read_csv(file.path(a, 'patient_edge_scores.csv.gz'), show_col_types = FALSE, progress = FALSE,
                col_types = cols_only(patient_id = col_character(), source_idx = col_integer(),
                                      target_idx = col_integer(), median_gated_attention_change = col_double())) %>%
    filter(source_idx != target_idx) %>% mutate(u = pmin(source_idx, target_idx), v = pmax(source_idx, target_idx)) %>%
    group_by(patient_id, u, v) %>% summarise(change = mean(median_gated_attention_change), .groups = 'drop') %>%
    left_join(meta, by = 'patient_id') %>%
    group_by(u, v) %>% summarise(edge_score = abs(rho(change, exposure)), .groups = 'drop') %>%
    left_join(rename(truth, u = node_idx, su = set_membership), by = 'u') %>%
    left_join(rename(truth, v = node_idx, sv = set_membership), by = 'v') %>%
    mutate(seed = id, edge_class = if_else(su == sv, su, 'Mixed'))
  edge_scores[[i]] <- e
  # Weighted degree of each node over the top 10% edges (within-class edges only).
  top <- e %>% filter(is.finite(edge_score), edge_score > 0) %>% arrange(desc(edge_score), u, v) %>%
    slice_head(n = ceiling(0.10 * nrow(e))) %>% filter(edge_class != 'Mixed')
  w <- bind_rows(transmute(top, node_idx = u, edge_score), transmute(top, node_idx = v, edge_score)) %>%
    group_by(node_idx) %>% summarise(weighted_degree = sum(edge_score), .groups = 'drop')
  degree_rows[[i]] <- truth %>% left_join(w, by = 'node_idx') %>%
    mutate(weighted_degree = replace_na(weighted_degree, 0), seed = id)
}
edge_scores <- bind_rows(edge_scores) %>% filter(edge_class != 'Mixed', is.finite(edge_score))

# A. Density of edge scores by class: estimated per dataset (unit area, reflection at 0 and 1)
#    and averaged equally across datasets.
bw <- 0.05; xs <- seq(0, 1, length.out = 401)
dens <- edge_scores %>% group_by(seed, edge_class) %>% group_modify(~{
  v <- .x$edge_score
  y <- vapply(xs, function(x) mean(dnorm((x - v) / bw) + dnorm((x + v) / bw) + dnorm((x - (2 - v)) / bw)) / bw, numeric(1))
  tibble(edge_score = xs, density = y / sum((head(y, -1) + tail(y, -1)) / 2 * diff(xs)))
}) %>% group_by(edge_class, edge_score) %>% summarise(density = mean(density), .groups = 'drop') %>%
  mutate(edge_class = factor(edge_class, levels = c('A', 'C', 'B')))
write_csv(dens, file.path(fig_dir, 'Supplementary_Figure_05A_data.csv'))
edge_labels <- c(A = 'Set A (A-A)', C = 'Set C (C-C)', B = 'Set B (B-B)')
s5a <- ggplot(dens, aes(edge_score, density, fill = edge_class, color = edge_class)) +
  geom_area(position = 'identity', alpha = 0.28, linewidth = 0.65) +
  scale_fill_manual(values = set_colors, labels = edge_labels) + scale_color_manual(values = set_colors, labels = edge_labels) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Exposure association of edge attention changes', x = 'Edge score', y = 'Density',
       fill = 'Edge class', color = 'Edge class')
save_plot(s5a, 'Supplementary_Figure_05A', 8, 5)

# B. Mean weighted degree per class: averaged within each dataset; bars = mean, error bars = SD across datasets.
seed_degree <- bind_rows(degree_rows) %>% group_by(seed, set_membership) %>%
  summarise(weighted_degree = mean(weighted_degree), .groups = 'drop') %>%
  mutate(set = factor(set_membership, levels = c('A', 'C', 'B')))
write_csv(seed_degree, file.path(fig_dir, 'Supplementary_Figure_05B_data.csv'))
deg_summary <- seed_degree %>% group_by(set) %>% summarise(mean = mean(weighted_degree), sd = sd(weighted_degree), .groups = 'drop')
s5b <- ggplot(deg_summary, aes(set, mean, fill = set)) + geom_col(width = 0.68, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.18, na.rm = TRUE) +
  geom_point(data = seed_degree, aes(y = weighted_degree), position = position_jitter(width = 0.055, height = 0), size = 1.3, alpha = 0.65) +
  geom_text(aes(y = mean + coalesce(sd, 0), label = round(mean, 2)), vjust = -0.45) +
  scale_fill_manual(values = set_colors) + scale_x_discrete(labels = set_labels) +
  labs(title = 'Weighted degree among top 10% exposure-associated edges', x = NULL, y = 'Average weighted degree') +
  theme(legend.position = 'none')
save_plot(s5b, 'Supplementary_Figure_05B', 7, 5.5)

message('Simulation figures written to ', normalizePath(fig_dir))
