# =============================================================================
# Node-score heatmaps across training lengths (genes and CpGs separately)
# =============================================================================
# Output: Heatmap_genes_<VALUE>_by_epoch.pdf and Heatmap_CpGs_<VALUE>_by_epoch.pdf, one page per
# analysis run (e.g. 400, 600, 800 epochs), plus PNG pages and CSV tables.
#
# Run from Terminal:
#   Rscript heatmap_epochs.R INPUT_DATA OUTPUT_DIR ANALYSIS_DIR_1 ANALYSIS_DIR_2 [ANALYSIS_DIR_3 ...]
# Example:
#   Rscript heatmap_epochs.R "$C" "$V/heatmaps_by_epoch" \
#     "$V/TCGA_clean59_e400/analysis" "$V/TCGA_clean59_e600/analysis" "$V/TCGA_clean59_e800/analysis"
#
# Design (keeps pages comparable and avoids picking the best-looking run):
# - The SAME nodes, in the SAME order, appear on every page. They are the top nodes
#   of the reference run (default: the 600-epoch run; set REFERENCE_EPOCHS=... to change).
# - Columns = patients ordered by recorded exposure; top strip shows exposure.
# - Cells = patient-specific node score (median |contribution| across initializations),
#   z-scored within each node so exposure patterns are visible.
# - Right bar = each node's Spearman correlation with recorded exposure (descriptive).
# - Each page also reports alignment with shuffled labels: permutation runs scored
#   against THEIR OWN shuffled exposure. Real-minus-shuffled is the alignment beyond
#   what fitting alone produces. Longer training raises both; compare the gap.
# No p-values; in-sample, descriptive. Model scores are read, never changed.
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
top_n_genes <- 40        # gene rows per page
top_n_cpgs <- 30         # CpG rows per page (all CpGs if fewer)
reference_epochs <- as.integer(Sys.getenv('REFERENCE_EPOCHS', '600'))
# VALUE=priority (default): |contribution|, the score analyze.py ranks nodes by.
# VALUE=signed: signed contribution (direction of push on the exposure score).
# The fitted exposure score is the SUM of signed contributions, so exposure
# grading is expected in signed values; |contribution| loses the direction.
value_type <- Sys.getenv('VALUE', 'priority')
if (!value_type %in% c('priority', 'signed')) stop('VALUE must be priority or signed')
npz_key <- if (value_type == 'signed') 'signed' else 'score'

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('Usage: Rscript heatmap_epochs.R INPUT_DATA OUTPUT_DIR ANALYSIS_DIR_1 [ANALYSIS_DIR_2 ...]')
data_dir <- args[1]; output_dir <- args[2]; analysis_dirs <- args[-(1:2)]
for (d in c(data_dir, analysis_dirs)) if (!dir.exists(d)) stop('Folder not found: ', d)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) + theme(plot.title = element_text(face = 'bold'),
                                           panel.grid = element_blank()))
read_table <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)

# =============================================================================
# Helpers
# =============================================================================
# Reads one array from numpy .npz (analyze.py attributions.npz) without Python.
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
# Spearman of every column of m with y (NA for constant columns).
column_spearman <- function(m, y) {
  ry <- rank(y)
  rm <- apply(m, 2, rank)
  suppressWarnings(as.vector(cor(rm, ry)))
}

# =============================================================================
# Recorded exposure for the fitted cohort
# =============================================================================
first_config <- fromJSON(file.path(analysis_dirs[1], 'config.json'))
exposure_column <- first_config$exposure_col
meta <- read_table(file.path(data_dir, 'metadata_samples.csv')) %>%
  transmute(patient_id = as.character(patient_id), exposure = .data[[exposure_column]]) %>% distinct()
if (anyDuplicated(meta$patient_id)) stop('Different exposure values for the same patient.')

# =============================================================================
# Load each run
# =============================================================================
runs <- list()
for (a in analysis_dirs) {
  config <- fromJSON(file.path(a, 'config.json'))
  label <- paste0(config$epochs, ' epochs')
  scores <- read_table(file.path(a, 'patient_node_scores.csv')) %>%
    mutate(patient_id = as.character(patient_id)) %>%
    select(patient_id, node_idx, node_label, median_priority)
  cohort <- read_table(file.path(a, 'cohort_summary.csv')) %>% select(node_idx, node_label, median_patient_rank)
  missing <- setdiff(unique(scores$patient_id), meta$patient_id)
  if (length(missing)) stop('Patients in ', a, ' missing from INPUT_DATA metadata: ', paste(head(missing), collapse = ', '))
  run_table <- read_table(file.path(a, 'runs.csv')) %>%
    filter(kind %in% c('initialization', 'permutation')) %>%
    mutate(folder = file.path(a, paste0(kind, '_', .data[['repeat']])))
  if (value_type == 'signed') {
    init <- run_table %>% filter(kind == 'initialization')
    ids <- as.character(read_table(file.path(init$folder[1], 'fit_membership.csv'))$patient_id)
    stack <- simplify2array(lapply(init$folder, function(f) read_npy_from_npz(file.path(f, 'attributions.npz'), 'signed')))
    med <- apply(stack, c(1, 2), median)
    labels_by_idx <- scores %>% distinct(node_idx, node_label) %>% arrange(node_idx)
    scores <- tibble(patient_id = rep(ids, times = ncol(med)), node_idx = rep(seq_len(ncol(med)) - 1L, each = nrow(med)),
                     median_priority = as.vector(med)) %>%
      left_join(labels_by_idx, by = 'node_idx') %>% select(patient_id, node_idx, node_label, median_priority)
  }
  runs[[label]] <- list(dir = a, epochs = config$epochs, label = label, scores = scores,
                        cohort = cohort, run_table = run_table)
}
if (anyDuplicated(names(runs))) stop('Two analysis folders have the same epoch count.')
epoch_values <- sapply(runs, `[[`, 'epochs')
reference <- if (reference_epochs %in% epoch_values) names(runs)[epoch_values == reference_epochs] else names(runs)[1]
message('Reference run for node selection and order: ', reference)

# Per-run alignment: per-node Spearman of patient scores with the supplied exposure,
# for each initialization (real labels) and permutation (own shuffled labels) run.
alignment_by_run <- function(run) {
  out <- list()
  for (j in seq_len(nrow(run$run_table))) {
    folder <- run$run_table$folder[j]
    npz <- file.path(folder, 'attributions.npz')
    if (!file.exists(npz)) next
    membership <- read_table(file.path(folder, 'fit_membership.csv')) %>% mutate(patient_id = as.character(patient_id))
    score <- read_npy_from_npz(npz, npz_key)
    out[[length(out) + 1]] <- tibble(kind = run$run_table$kind[j], repeat_id = run$run_table[['repeat']][j],
      node_idx = seq_len(ncol(score)) - 1L, rho = column_spearman(score, membership$supplied_exposure))
  }
  bind_rows(out)
}
for (nm in names(runs)) {
  message('Computing real vs shuffled alignment for ', nm)
  runs[[nm]]$alignment <- alignment_by_run(runs[[nm]])
}

# =============================================================================
# Heatmap pages
# =============================================================================
make_page <- function(run, nodes_shown, type_label, all_type_idx) {
  pat <- meta %>% filter(patient_id %in% run$scores$patient_id) %>% arrange(exposure, patient_id)
  dat <- run$scores %>% filter(node_idx %in% nodes_shown$node_idx) %>%
    left_join(pat, by = 'patient_id') %>% group_by(node_idx) %>%
    mutate(z = if (sd(median_priority) > 0) (median_priority - mean(median_priority)) / sd(median_priority) else 0) %>%
    ungroup() %>%
    mutate(z = pmax(pmin(z, 3), -3),
           patient_id = factor(patient_id, levels = pat$patient_id),
           display = factor(nodes_shown$display[match(node_idx, nodes_shown$node_idx)], levels = rev(nodes_shown$display)))
  node_rho <- dat %>% group_by(node_idx, display) %>%
    summarise(rho = suppressWarnings(cor(rank(median_priority), rank(exposure))), .groups = 'drop')

  # Alignment summaries for the shown nodes and for all nodes of this type.
  al <- run$alignment
  summarise_alignment <- function(idx) {
    al %>% filter(node_idx %in% idx) %>% group_by(kind, repeat_id) %>%
      summarise(m = median(abs(rho), na.rm = TRUE), .groups = 'drop') %>%
      group_by(kind) %>% summarise(value = median(m), n_runs = n(), .groups = 'drop')
  }
  shown_al <- summarise_alignment(nodes_shown$node_idx); all_al <- summarise_alignment(all_type_idx)
  get <- function(t, k) { v <- t$value[t$kind == k]; if (length(v)) sprintf('%.3f', v) else 'NA' }
  n_init <- sum(al$kind == 'initialization' & al$node_idx == al$node_idx[1])
  n_perm <- sum(al$kind == 'permutation' & al$node_idx == al$node_idx[1])
  subtitle <- paste0('Median |Spearman| of node score with exposure (per run, then median over runs)\n',
    'Shown nodes: real labels ', get(shown_al, 'initialization'), ' vs shuffled (own labels) ', get(shown_al, 'permutation'),
    '   |   All ', type_label, ': real ', get(all_al, 'initialization'), ' vs shuffled ', get(all_al, 'permutation'),
    '\n', n_init, ' initialization and ', n_perm, ' permutation runs; ', nrow(pat), ' patients')

  strip <- ggplot(pat %>% mutate(patient_id = factor(patient_id, levels = pat$patient_id)),
                  aes(patient_id, 'Exposure', fill = exposure)) + geom_tile() +
    scale_fill_gradient(low = 'white', high = '#3B3B3B', name = paste0('Recorded\n', exposure_column)) +
    labs(x = NULL, y = NULL, title = paste0(type_label, ': ', run$label, if (value_type == 'signed') ' (signed contribution)' else ' (node priority |contribution|)'),
         subtitle = subtitle) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(size = 8.5))
  heat <- ggplot(dat, aes(patient_id, display, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = '#25259A', mid = 'white', high = '#E52727', midpoint = 0,
                         limits = c(-3, 3), name = if (value_type == 'signed') 'Signed contribution\n(z within node)' else 'Node score\n(z within node)') +
    labs(x = 'Patients ordered by recorded exposure (low to high)', y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7))
  bar <- ggplot(node_rho, aes(rho, display, fill = rho > 0)) + geom_col(width = 0.75) +
    geom_vline(xintercept = 0, color = 'grey40') +
    scale_fill_manual(values = c(`TRUE` = '#E52727', `FALSE` = '#25259A'), guide = 'none') +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    labs(x = 'Spearman\nwith exposure', y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  page <- (strip + plot_spacer() + plot_layout(widths = c(1, 0.18))) /
    (heat + bar + plot_layout(widths = c(1, 0.18))) +
    plot_layout(heights = c(0.07, 1), guides = 'collect') +
    plot_annotation(caption = paste0('Same ', type_label, ' and order on every page: top ', nrow(nodes_shown),
      ' by cohort median rank in the ', reference, ' run, ordered by their Spearman there. Descriptive, in-sample; no p-values.'))
  list(page = page, rho = node_rho %>% mutate(run = run$label),
       alignment = bind_rows(shown_al %>% mutate(nodes = 'shown'), all_al %>% mutate(nodes = paste('all', type_label))) %>%
         mutate(run = run$label, type = type_label))
}

summary_rows <- list()
for (type in c('GENE', 'CpG')) {
  type_label <- if (type == 'GENE') 'genes' else 'CpGs'
  top_n <- if (type == 'GENE') top_n_genes else top_n_cpgs
  ref <- runs[[reference]]
  type_nodes <- ref$cohort %>% filter(startsWith(node_label, paste0(type, ':')))
  if (!nrow(type_nodes)) { message('No ', type_label, ' found; skipping.'); next }
  selected <- type_nodes %>% arrange(median_patient_rank, node_idx) %>% slice_head(n = top_n)
  # Order rows by the reference run's node-exposure Spearman so patterns read top to bottom.
  ref_rho <- ref$scores %>% filter(node_idx %in% selected$node_idx) %>% left_join(meta, by = 'patient_id') %>%
    group_by(node_idx) %>% summarise(rho = suppressWarnings(cor(rank(median_priority), rank(exposure))), .groups = 'drop')
  selected <- selected %>% left_join(ref_rho, by = 'node_idx') %>% arrange(desc(rho), node_idx) %>%
    mutate(display = make.unique(sub('^(GENE|CpG):', '', node_label)))
  pdf_file <- file.path(output_dir, paste0('Heatmap_', type_label, '_', value_type, '_by_epoch.pdf'))
  height <- max(6, 2.6 + 0.16 * nrow(selected))
  pdf(pdf_file, width = 12, height = height, onefile = TRUE)
  rho_all <- list()
  for (nm in names(runs)[order(epoch_values)]) {
    res <- make_page(runs[[nm]], selected, type_label, type_nodes$node_idx)
    print(res$page)
    ggsave(file.path(output_dir, paste0('Heatmap_', type_label, '_', value_type, '_', runs[[nm]]$epochs, 'epochs.png')),
           res$page, width = 12, height = height, dpi = 200)
    rho_all[[nm]] <- res$rho; summary_rows[[length(summary_rows) + 1]] <- res$alignment
  }
  dev.off()
  bind_rows(rho_all) %>% select(run, node_idx, node = display, spearman_with_exposure = rho) %>%
    pivot_wider(names_from = run, values_from = spearman_with_exposure) %>%
    write_csv(file.path(output_dir, paste0('Heatmap_', type_label, '_', value_type, '_node_spearman_by_epoch.csv')))
  message('Saved ', pdf_file)
}
alignment <- bind_rows(summary_rows) %>%
  select(type, run, nodes, kind, median_abs_spearman = value, n_runs) %>%
  pivot_wider(names_from = kind, values_from = c(median_abs_spearman, n_runs)) %>%
  mutate(real_minus_shuffled = median_abs_spearman_initialization - median_abs_spearman_permutation)
write_csv(alignment, file.path(output_dir, paste0('Heatmap_alignment_real_vs_shuffled_', value_type, '.csv')))
print(as.data.frame(alignment))
writeLines(c(
  paste('Reference run for node selection/order:', reference),
  paste('Cell values:', if (value_type == 'signed') 'signed contribution (median over initializations)' else 'node priority = median |contribution| over initializations'),
  'Rows are identical across pages so runs can be compared; they are not re-selected per run.',
  'Alignment = median |Spearman| between patient node score and exposure. Real runs use recorded',
  'exposure; permutation runs use their own shuffled exposure. Longer training raises both;',
  'real_minus_shuffled is the alignment beyond fitting. Descriptive, in-sample, no p-values.',
  'Use these pages to check robustness across training lengths, not to choose the run after viewing.'
), file.path(output_dir, paste0('Heatmap_notes_', value_type, '.txt')))
