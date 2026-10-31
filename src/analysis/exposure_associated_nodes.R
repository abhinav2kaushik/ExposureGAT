# =============================================================================
# Exposure-associated node score and heatmaps (genes and CpGs separately)
# =============================================================================
# Post-hoc; reads ExposureGAT model outputs only; model scores are unchanged.
#
# Node priority (model output) = how much a node contributes to the fitted exposure
# score. A node can matter to the fit without tracking exposure. This script adds:
#
#   exposure alignment  rho(node) = Spearman between the node's SIGNED patient
#                       contribution and recorded exposure (median over the
#                       initialization runs).
#   null percentile     fraction of shuffled-label |rho| values (permutation runs,
#                       each against ITS OWN shuffled exposure, pooled over runs and
#                       over nodes of the same type) that |rho(node)| exceeds.
#                       This is what fitting alone produces, so it calibrates rho.
#   exposure-associated node score = normalized node priority x null percentile
#                       (cohort level; patient level uses the patient priority).
#   exposure-associated node = null percentile >= THRESHOLD (default 0.99).
#
# Checks written alongside:
#   - observed number of passing nodes vs the number expected by chance, and the
#     same count in every real and every shuffled run (same rule, same null);
#   - overlap with a plain input-data correlation ranking (does the GNN add anything?).
# Descriptive, in-sample; null percentiles are not p-values.
#
# Run from Terminal (one or more analysis folders; each gets its own page):
#   Rscript exposure_associated_nodes.R INPUT_DATA OUTPUT_DIR ANALYSIS_DIR [ANALYSIS_DIR ...]
# Example:
#   Rscript exposure_associated_nodes.R "$V/GNN_Input_Data_clean59" "$V/exposure_associated" \
#     "$V/TCGA_clean59_e600/analysis"
# Options (environment variables): THRESHOLD=0.99 TOP_GENES=40 TOP_CPGS=30
# install.packages(c('dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite'))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)

# =============================================================================
# Settings
# =============================================================================
threshold <- as.numeric(Sys.getenv('THRESHOLD', '0.99'))
top_n <- c(GENE = as.integer(Sys.getenv('TOP_GENES', '40')), CpG = as.integer(Sys.getenv('TOP_CPGS', '30')))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('Usage: Rscript exposure_associated_nodes.R INPUT_DATA OUTPUT_DIR ANALYSIS_DIR [ANALYSIS_DIR ...]')
data_dir <- args[1]; output_dir <- args[2]; analysis_dirs <- args[-(1:2)]
for (d in c(data_dir, analysis_dirs)) if (!dir.exists(d)) stop('Folder not found: ', d)
if (!(threshold > 0 && threshold < 1)) stop('THRESHOLD must be between 0 and 1')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) + theme(plot.title = element_text(face = 'bold')))
read_table <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)
type_names <- c(GENE = 'genes', CpG = 'CpGs')

# =============================================================================
# Helpers
# =============================================================================
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
column_spearman <- function(m, y) suppressWarnings(as.vector(cor(apply(m, 2, rank), rank(y))))
# Fraction of null values strictly below each observed value (ties count half).
null_percentile <- function(observed, null) {
  null <- sort(null[is.finite(null)])
  if (!length(null)) return(rep(NA_real_, length(observed)))
  below <- findInterval(observed, null, left.open = TRUE)
  at_or_below <- findInterval(observed, null)
  ifelse(is.finite(observed), (below + at_or_below) / 2 / length(null), NA_real_)
}

# =============================================================================
# Recorded exposure and paired input changes (for the plain-correlation check)
# =============================================================================
first_config <- fromJSON(file.path(analysis_dirs[1], 'config.json'))
exposure_column <- first_config$exposure_col
sample_meta <- read_table(file.path(data_dir, 'metadata_samples.csv')) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id))
meta <- sample_meta %>% transmute(patient_id, exposure = .data[[exposure_column]]) %>% distinct()
if (anyDuplicated(meta$patient_id)) stop('Different exposure values for the same patient.')
pairs <- inner_join(sample_meta %>% filter(grepl('NT', tissue_type, ignore.case = TRUE)) %>% select(patient_id, normal_id = sample_id),
                    sample_meta %>% filter(grepl('TP', tissue_type, ignore.case = TRUE)) %>% select(patient_id, tumor_id = sample_id),
                    by = 'patient_id') %>% distinct(patient_id, .keep_all = TRUE) %>% left_join(meta, by = 'patient_id')
input_delta_rho <- function(file, prefix) {
  x <- read_table(file.path(data_dir, file)); names(x)[1] <- 'sample_id'
  x$sample_id <- as.character(x$sample_id)
  keep <- pairs %>% filter(normal_id %in% x$sample_id, tumor_id %in% x$sample_id)
  xm <- as.matrix(x[, -1])
  delta <- xm[match(keep$tumor_id, x$sample_id), , drop = FALSE] - xm[match(keep$normal_id, x$sample_id), , drop = FALSE]
  tibble(node_label = paste0(prefix, colnames(xm)), rho = column_spearman(delta, keep$exposure))
}
message('Computing plain input-data correlations (paired tumor - normal)')
plain <- bind_rows(
  input_delta_rho('X_rna_gene_features.csv', 'GENE:') %>% mutate(modality = 'rna'),
  input_delta_rho('X_meth_gene_features.csv', 'GENE:') %>% mutate(modality = 'meth'),
  input_delta_rho('X_cnv_gene_features.csv', 'GENE:') %>% mutate(modality = 'cnv'),
  input_delta_rho('X_cpg_features.csv', 'CpG:') %>% mutate(modality = 'beta')) %>%
  group_by(node_label) %>% slice_max(abs(rho), n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(node_label, input_delta_rho = rho, input_delta_modality = modality)

# =============================================================================
# One analysis run
# =============================================================================
analyse_run <- function(a) {
  config <- fromJSON(file.path(a, 'config.json'))
  label <- paste0(config$epochs, ' epochs')
  run_table <- read_table(file.path(a, 'runs.csv')) %>% filter(kind %in% c('initialization', 'permutation')) %>%
    mutate(folder = file.path(a, paste0(kind, '_', .data[['repeat']])))
  if (sum(run_table$kind == 'initialization') < 1 || sum(run_table$kind == 'permutation') < 2)
    stop('Need initialization and >=2 permutation runs in ', a)
  cohort <- read_table(file.path(a, 'cohort_summary.csv')) %>% arrange(node_idx)
  n_nodes <- nrow(cohort)
  node_type <- ifelse(startsWith(cohort$node_label, 'CpG:'), 'CpG', 'GENE')
  ids <- NULL; rho <- list(); signed_init <- list()
  for (j in seq_len(nrow(run_table))) {
    f <- run_table$folder[j]
    m <- read_table(file.path(f, 'fit_membership.csv')) %>% mutate(patient_id = as.character(patient_id))
    if (is.null(ids)) ids <- m$patient_id else if (!identical(ids, m$patient_id)) stop('Patient order differs in ', f)
    signed <- read_npy_from_npz(file.path(f, 'attributions.npz'), 'signed')
    if (run_table$kind[j] == 'initialization') {
      recorded <- meta$exposure[match(ids, meta$patient_id)]
      if (any(abs(recorded - m$supplied_exposure) > 1e-6)) stop('INPUT_DATA exposure differs from the fitted exposure; use the folder the run was fitted on.')
      signed_init[[length(signed_init) + 1]] <- signed
    }
    rho[[j]] <- column_spearman(signed, m$supplied_exposure)
  }
  rho <- do.call(rbind, rho)                                   # runs x nodes
  is_init <- run_table$kind == 'initialization'
  real_rho <- apply(rho[is_init, , drop = FALSE], 2, median)

  # Pooled shuffled-label null per node type; percentiles for real nodes.
  pct <- rep(NA_real_, n_nodes); counts <- list()
  for (t in c('GENE', 'CpG')) {
    idx <- which(node_type == t); if (!length(idx)) next
    null_all <- abs(rho[!is_init, idx, drop = FALSE])
    pct[idx] <- null_percentile(abs(real_rho[idx]), as.vector(null_all))
    # Same rule applied to every individual run (shuffled runs: leave-one-out null).
    perm_rows <- which(!is_init)
    for (j in seq_len(nrow(rho))) {
      null_j <- if (is_init[j]) null_all else abs(rho[setdiff(perm_rows, j), idx, drop = FALSE])
      p_j <- null_percentile(abs(rho[j, idx]), as.vector(null_j))
      counts[[length(counts) + 1]] <- tibble(type = type_names[[t]], kind = run_table$kind[j],
        repeat_id = run_table[['repeat']][j], n_passing = sum(p_j >= threshold, na.rm = TRUE), n_nodes = length(idx))
    }
  }
  counts <- bind_rows(counts) %>% mutate(run = label)

  nodes <- cohort %>% transmute(node_idx, node_label, type = type_names[node_type], median_patient_rank,
      normalized_priority = (n_nodes - median_patient_rank) / max(n_nodes - 1, 1)) %>%
    mutate(exposure_rho = real_rho,
           exposure_rho_min = apply(rho[is_init, , drop = FALSE], 2, min),
           exposure_rho_max = apply(rho[is_init, , drop = FALSE], 2, max),
           null_percentile = pct,
           exposure_associated_score = normalized_priority * null_percentile,
           exposure_associated = null_percentile >= threshold) %>%
    left_join(plain, by = 'node_label') %>% mutate(run = label)
  summary <- nodes %>% group_by(run, type) %>% summarise(n_nodes = n(),
      observed_exposure_associated = sum(exposure_associated, na.rm = TRUE),
      expected_by_chance = round(n_nodes * (1 - threshold), 1), .groups = 'drop') %>%
    left_join(counts %>% group_by(type, kind) %>% summarise(median_count = median(n_passing), .groups = 'drop') %>%
                pivot_wider(names_from = kind, values_from = median_count, names_prefix = 'median_count_per_run_'),
              by = 'type') %>% mutate(threshold = threshold)

  # Patient-level exposure-associated score for every node.
  signed_med <- apply(simplify2array(signed_init), c(1, 2), median)
  patient <- read_table(file.path(a, 'patient_node_scores.csv')) %>% mutate(patient_id = as.character(patient_id)) %>%
    select(patient_id, node_idx, median_priority, median_rank) %>%
    mutate(patient_normalized_priority = (n_nodes - median_rank) / max(n_nodes - 1, 1),
           median_signed_contribution = signed_med[cbind(match(patient_id, ids), node_idx + 1L)]) %>%
    left_join(nodes %>% select(node_idx, node_label, type, null_percentile, exposure_associated), by = 'node_idx') %>%
    mutate(patient_exposure_associated_score = patient_normalized_priority * null_percentile) %>%
    left_join(meta, by = 'patient_id')
  list(label = label, epochs = config$epochs, nodes = nodes, counts = counts, summary = summary, patient = patient,
       n_init = sum(is_init), n_perm = sum(!is_init))
}

results <- lapply(analysis_dirs, function(a) { message('Analysing ', a); analyse_run(a) })
names(results) <- sapply(results, `[[`, 'label')
results <- results[order(sapply(results, `[[`, 'epochs'))]

# =============================================================================
# Tables
# =============================================================================
for (r in results) {
  tag <- paste0(r$epochs, 'epochs')
  write_csv(r$nodes %>% arrange(type, desc(exposure_associated_score)), file.path(output_dir, paste0('Exposure_associated_node_scores_', tag, '.csv')))
  write_csv(r$patient %>% filter(exposure_associated) %>% arrange(type, node_idx, patient_id),
            file.path(output_dir, paste0('Exposure_associated_patient_scores_', tag, '.csv')))
}
summary_all <- bind_rows(lapply(results, `[[`, 'summary'))
write_csv(summary_all, file.path(output_dir, 'Exposure_associated_summary.csv'))
print(as.data.frame(summary_all))

# Overlap with the original priority top list and with plain input correlation.
overlap <- bind_rows(lapply(results, function(r) r$nodes %>% group_by(run, type) %>% group_modify(function(d, key) {
  sel <- d %>% filter(exposure_associated) %>% pull(node_idx); k <- length(sel)
  if (!k) return(tibble(n_selected = 0L))
  top_priority <- d %>% arrange(median_patient_rank) %>% slice_head(n = k) %>% pull(node_idx)
  top_plain <- d %>% arrange(desc(abs(input_delta_rho))) %>% slice_head(n = k) %>% pull(node_idx)
  tibble(n_selected = k, shared_with_priority_top_k = length(intersect(sel, top_priority)),
         shared_with_plain_input_correlation_top_k = length(intersect(sel, top_plain)),
         expected_shared_by_chance = round(k * k / nrow(d), 1))
}) %>% ungroup()))
write_csv(overlap, file.path(output_dir, 'Exposure_associated_overlap_checks.csv'))

# =============================================================================
# Figure: observed vs chance, and the same rule in real vs shuffled runs
# =============================================================================
counts_all <- bind_rows(lapply(results, `[[`, 'counts')) %>%
  mutate(kind = factor(recode(kind, initialization = 'Real labels', permutation = 'Shuffled labels'),
                       levels = c('Real labels', 'Shuffled labels')),
         run = factor(run, levels = names(results)))
p_counts <- ggplot(counts_all, aes(kind, n_passing, color = kind)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = 'grey45') +
  geom_jitter(width = 0.1, height = 0, size = 1.8, alpha = 0.8) +
  geom_hline(data = summary_all %>% mutate(run = factor(run, levels = names(results))),
             aes(yintercept = expected_by_chance), linetype = 'dashed') +
  facet_grid(type ~ run, scales = 'free_y') +
  scale_color_manual(values = c('#2886B5', '#A14792')) +
  labs(title = 'Nodes passing the exposure-alignment rule, per run',
       subtitle = paste0('Rule: |Spearman(signed contribution, exposure)| above ', threshold * 100,
                         '% of shuffled-label values; dashed = expected by chance'),
       x = NULL, y = 'Nodes passing',
       caption = 'Each dot is one fit. Shuffled runs are scored against their own shuffled exposure.') +
  theme(legend.position = 'none')
ggsave(file.path(output_dir, 'Exposure_associated_counts_real_vs_shuffled.pdf'), p_counts, width = 4 + 3 * length(results), height = 7)
ggsave(file.path(output_dir, 'Exposure_associated_counts_real_vs_shuffled.png'), p_counts, width = 4 + 3 * length(results), height = 7, dpi = 200)

# =============================================================================
# Heatmaps: exposure-associated nodes only (one PDF per node type, page per run)
# =============================================================================
make_page <- function(r, t) {
  sel <- r$nodes %>% filter(type == type_names[[t]], exposure_associated) %>%
    arrange(desc(exposure_associated_score), node_idx) %>% slice_head(n = top_n[[t]])
  n_pass <- sum(r$nodes$type == type_names[[t]] & r$nodes$exposure_associated, na.rm = TRUE)
  expected <- r$summary$expected_by_chance[r$summary$type == type_names[[t]]]
  head_text <- paste0(type_names[[t]], ': exposure-associated nodes, ', r$label)
  sub_text <- paste0(n_pass, ' passing (expected by chance ', expected, '); showing top ', nrow(sel),
                     ' by exposure-associated score. ', r$n_init, ' initializations, ', r$n_perm,
                     ' permutations; ', n_distinct(r$patient$patient_id), ' patients.')
  if (!nrow(sel)) {
    return(ggplot() + annotate('text', x = 0, y = 0, label = 'No nodes passed the exposure-alignment rule.', size = 6) +
             theme_void() + labs(title = head_text, subtitle = sub_text))
  }
  sel <- sel %>% arrange(desc(exposure_rho)) %>% mutate(display = make.unique(sub('^(GENE|CpG):', '', node_label)))
  pat <- meta %>% filter(patient_id %in% r$patient$patient_id) %>% arrange(exposure, patient_id)
  dat <- r$patient %>% filter(node_idx %in% sel$node_idx) %>% group_by(node_idx) %>%
    mutate(z = if (sd(median_signed_contribution) > 0)
             (median_signed_contribution - mean(median_signed_contribution)) / sd(median_signed_contribution) else 0) %>%
    ungroup() %>% mutate(z = pmax(pmin(z, 3), -3), patient_id = factor(patient_id, levels = pat$patient_id),
                         display = factor(sel$display[match(node_idx, sel$node_idx)], levels = rev(sel$display)))
  side <- sel %>% mutate(display = factor(display, levels = rev(sel$display)))
  strip <- ggplot(pat %>% mutate(patient_id = factor(patient_id, levels = pat$patient_id)), aes(patient_id, 'Exposure', fill = exposure)) +
    geom_tile() + scale_fill_gradient(low = 'white', high = '#3B3B3B', name = paste0('Recorded\n', exposure_column)) +
    labs(x = NULL, y = NULL, title = head_text, subtitle = sub_text) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          plot.subtitle = element_text(size = 8.5))
  heat <- ggplot(dat, aes(patient_id, display, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = '#25259A', mid = 'white', high = '#E52727', midpoint = 0, limits = c(-3, 3),
                         name = 'Signed contribution\n(z within node)') +
    labs(x = 'Patients ordered by recorded exposure (low to high)', y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          axis.text.y = element_text(size = 7))
  bar_rho <- ggplot(side, aes(exposure_rho, display, fill = exposure_rho > 0)) + geom_col(width = 0.75) +
    geom_errorbarh(aes(xmin = exposure_rho_min, xmax = exposure_rho_max), height = 0.3, color = 'grey30') +
    geom_vline(xintercept = 0, color = 'grey40') +
    scale_fill_manual(values = c(`TRUE` = '#E52727', `FALSE` = '#25259A'), guide = 'none') +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    labs(x = 'Spearman\nwith exposure', y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
  bar_score <- ggplot(side, aes(exposure_associated_score, display)) + geom_col(width = 0.75, fill = '#337A96') +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = 'Exposure-\nassociated score', y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
  (strip + plot_spacer() + plot_spacer() + plot_layout(widths = c(1, 0.16, 0.16))) /
    (heat + bar_rho + bar_score + plot_layout(widths = c(1, 0.16, 0.16))) +
    plot_layout(heights = c(0.07, 1), guides = 'collect') +
    plot_annotation(caption = paste0('Rows: nodes whose signed contribution tracks exposure beyond ', threshold * 100,
      '% of shuffled-label alignment, ordered by Spearman. Bars: median over initializations (whiskers = range).',
      '\nScore = node priority x null percentile. Descriptive, in-sample; not p-values.'))
}
for (t in c('GENE', 'CpG')) {
  if (!any(sapply(results, function(r) any(r$nodes$type == type_names[[t]])))) next
  n_max <- max(sapply(results, function(r) min(top_n[[t]], sum(r$nodes$type == type_names[[t]] & r$nodes$exposure_associated, na.rm = TRUE))))
  height <- max(6, 2.8 + 0.16 * n_max)
  pdf_file <- file.path(output_dir, paste0('Exposure_associated_', type_names[[t]], '_heatmap.pdf'))
  pdf(pdf_file, width = 13, height = height, onefile = TRUE)
  for (r in results) {
    page <- make_page(r, t); print(page)
    ggsave(file.path(output_dir, paste0('Exposure_associated_', type_names[[t]], '_heatmap_', r$epochs, 'epochs.png')),
           page, width = 13, height = height, dpi = 200)
  }
  dev.off(); message('Saved ', pdf_file)
}
writeLines(c(
  paste('Threshold (null percentile):', threshold),
  'exposure_rho: Spearman of signed node contribution with recorded exposure, median over initialization runs.',
  'null_percentile: share of pooled shuffled-label |rho| (permutation runs vs their own labels, same node type) below |exposure_rho|.',
  'exposure_associated_score = normalized cohort priority x null_percentile; patient score uses the patient priority.',
  'Compare observed_exposure_associated with expected_by_chance and with the shuffled-run counts before interpreting.',
  'If selected nodes mostly match the plain input-correlation top list, the GNN adds little for those nodes.',
  'Descriptive and in-sample. Null percentiles are not p-values and are not multiple-testing corrected.'
), file.path(output_dir, 'Exposure_associated_notes.txt'))
