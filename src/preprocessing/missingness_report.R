# =============================================================================
# Missingness report for the GNN inputs (values for the response and methods)
# =============================================================================
# Reconstructs, without rebuilding anything, how many values Step_01 replaced:
#   - methylation samples entirely missing  -> normal-tissue median
#   - CNV samples entirely missing          -> diploid
#   - individual missing CNV values         -> diploid
#   - gene nodes without promoter probes    -> within-sample median methylation
#   - individual missing CpG beta values    -> k-nearest-neighbour imputation
#   - exposure: zero / missing pack-years by smoking status
# Reported for all paired patients and for the analyzed cohort.
#
# Sample-level counts come from Session.RData using Step_01's matching rule
# (same tissue filter, 16-character barcode match). Gene-level methylation
# replacement is detected in the exported file: a gene with no promoter probe
# equals the within-sample median in every sample. Gene-level CNV replacement
# (genes absent from the CNV data) is detected as zero in every sample.
#
#   Rscript missingness_report.R SESSION_RDATA INPUT_DIR ANALYZED_INPUT_DIR OUTPUT_DIR
# Example:
#   Rscript missingness_report.R \
#     "$HOME/Library/CloudStorage/OneDrive-HarvardUniversity/GNN/TCGA/Lung/Session.RData" \
#     "$V/GNN_Input_Data_v2" "$V/GNN_Input_Data_v2_clean59" "$V/missingness"
# Packages: BiocManager::install('SummarizedExperiment')

suppressPackageStartupMessages(library(SummarizedExperiment))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('Usage: Rscript missingness_report.R SESSION_RDATA INPUT_DIR ANALYZED_INPUT_DIR OUTPUT_DIR')
session_file <- args[1]; input_dir <- args[2]; analyzed_dir <- args[3]; out_dir <- args[4]
for (p in c(input_dir, analyzed_dir)) if (!dir.exists(p)) stop('Folder not found: ', p)
if (!file.exists(session_file)) stop('Session file not found: ', session_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TARGET_TISSUES <- tolower(c('lower lobe, lung', 'lung, nos', 'main bronchus', 'middle lobe, lung', 'upper lobe, lung'))
read_matrix <- function(f) as.matrix(read.csv(file.path(input_dir, f), row.names = 1, check.names = FALSE))

# =============================================================================
# Exported inputs and cohorts
# =============================================================================
meta <- read.csv(file.path(input_dir, 'metadata_samples.csv'), check.names = FALSE)
sid <- intersect(c('sample_id', 'sample_barcode'), names(meta))[1]
meta$sample <- as.character(meta[[sid]]); meta$patient_id <- as.character(meta$patient_id)
analyzed_ids <- unique(as.character(read.csv(file.path(analyzed_dir, 'metadata_samples.csv'))$patient_id))
x_meth <- read_matrix('X_meth_gene_features.csv')
x_cnv  <- read_matrix('X_cnv_gene_features.csv')
x_cpg  <- read_matrix('X_cpg_features.csv')
samples <- rownames(x_meth)
cohorts <- list(all = unique(meta$patient_id), analyzed = analyzed_ids)
cat('Patients: all =', length(cohorts$all), ' analyzed =', length(cohorts$analyzed), '\n')

# =============================================================================
# Session.RData: Step_01 matching rule
# =============================================================================
env <- new.env(); load(session_file, envir = env)
for (o in c('data_meth', 'data_cnv')) if (!exists(o, envir = env)) stop(o, ' not found in ', session_file)
keep_tissue <- function(se) tolower(colData(se)$tissue_or_organ_of_origin) %in% TARGET_TISSUES
data_meth <- get('data_meth', envir = env); data_meth <- data_meth[, keep_tissue(data_meth)]
data_cnv  <- get('data_cnv',  envir = env); data_cnv  <- data_cnv[,  keep_tissue(data_cnv)]
rm(env); invisible(gc())
key16 <- substr(samples, 1, 16)
meth_idx <- match(key16, substr(colnames(data_meth), 1, 16))
cnv_idx  <- match(key16, substr(colnames(data_cnv), 1, 16))

# Individual missing CNV values among matched samples (all genes in the CNV data).
cnv_assay <- tryCatch(assay(data_cnv, 'copy_number'), error = function(e) assay(data_cnv))
# Individual missing CpG beta values for the exported CpG nodes, matched samples only.
cpg_ids <- intersect(colnames(x_cpg), rownames(data_meth))
beta_cpg <- assay(data_meth)[cpg_ids, , drop = FALSE]

# =============================================================================
# Gene-level replacements detected in exported files
# =============================================================================
sample_median <- apply(x_meth, 1, median)
meth_gene_replaced <- apply(abs(x_meth - sample_median) < 1e-9, 2, all)       # every sample = sample median
cnv_gene_replaced  <- apply(x_cnv == 0, 2, all)                               # zero in every sample

# =============================================================================
# Counts per cohort
# =============================================================================
pct <- function(n, d) if (d > 0) round(100 * n / d, 1) else NA
row <- function(item, replacement, unit, n, d) data.frame(item, replacement, unit, n, denominator = d, percent = pct(n, d))
count_cohort <- function(ids) {
  s <- meta$sample %in% samples & meta$patient_id %in% ids
  sm <- meta[s, ]; i <- match(sm$sample, samples)
  tumor <- sm$tissue_type == 'TP'
  m_ok <- !is.na(meth_idx[i]); c_ok <- !is.na(cnv_idx[i])
  cnv_vals <- if (any(c_ok)) cnv_assay[, cnv_idx[i][c_ok], drop = FALSE] else matrix(numeric(0))
  beta_vals <- if (any(m_ok)) beta_cpg[, meth_idx[i][m_ok], drop = FALSE] else matrix(numeric(0))
  pats <- sm[!duplicated(sm$patient_id), ]
  pk <- suppressWarnings(as.numeric(pats$smoking_intensity))
  status <- ifelse(is.na(pats$smoking_status) | pats$smoking_status == '', 'Not reported', as.character(pats$smoking_status))
  zero <- !is.na(pk) & pk == 0
  rbind(
    row('Patients', '', 'patients', nrow(pats), nrow(pats)),
    row('Methylation sample missing', 'normal-tissue median', 'samples', sum(!m_ok), nrow(sm)),
    row('  methylation missing, tumor samples', 'normal-tissue median', 'samples', sum(!m_ok & tumor), sum(tumor)),
    row('  methylation missing, normal samples', 'normal-tissue median', 'samples', sum(!m_ok & !tumor), sum(!tumor)),
    row('CNV sample missing', 'diploid (2)', 'samples', sum(!c_ok), nrow(sm)),
    row('  CNV missing, tumor samples', 'diploid (2)', 'samples', sum(!c_ok & tumor), sum(tumor)),
    row('  CNV missing, normal samples', 'diploid (2)', 'samples', sum(!c_ok & !tumor), sum(!tumor)),
    row('Individual CNV values missing (matched samples, all CNV genes)', 'diploid (2)', 'values',
        sum(is.na(cnv_vals)), length(cnv_vals)),
    row('Gene nodes without promoter methylation probes', 'within-sample median', 'gene nodes',
        sum(meth_gene_replaced), ncol(x_meth)),
    row('Gene nodes without CNV data (zero in all samples)', 'diploid (0 after centering)', 'gene nodes',
        sum(cnv_gene_replaced), ncol(x_cnv)),
    row('Individual CpG beta values missing (matched samples)', 'k-nearest-neighbour', 'values',
        sum(is.na(beta_vals)), length(beta_vals)),
    row('Exposure: zero pack-years recorded', 'exclusion unless lifelong non-smoker', 'patients', sum(zero), nrow(pats)),
    do.call(rbind, lapply(sort(unique(status[zero])), function(st)
      row(paste0('  zero pack-years, status: ', st), '', 'patients', sum(zero & status == st), sum(zero))))
  )
}
all_tab <- count_cohort(cohorts$all); an_tab <- count_cohort(cohorts$analyzed)
stopifnot(!anyDuplicated(all_tab$item), !anyDuplicated(an_tab$item))
out <- merge(cbind(order = seq_len(nrow(all_tab)), all_tab),
             an_tab[, c('item', 'n', 'denominator', 'percent')],
             by = 'item', all.x = TRUE, suffixes = c('_all', '_analyzed'), sort = FALSE)
out <- out[order(out$order), setdiff(names(out), 'order')]
write.csv(out, file.path(out_dir, 'Missingness_table.csv'), row.names = FALSE, na = '')
print(out, row.names = FALSE)
writeLines(c(
  'Missingness report',
  'Sample-level counts reproduce Step_01 matching (tissue filter + 16-character barcode).',
  'Gene-level methylation replacement: gene equal to the within-sample median in every sample.',
  'Gene-level CNV replacement: gene equal to zero (diploid after centering) in every sample;',
  '  a gene that is truly diploid in every tumor would also be counted, so this is an upper bound.',
  'Individual CNV missing values are counted over all genes in the CNV data, not only gene nodes.',
  'In the analyzed cohort, zero pack-years remain only for confirmed lifelong non-smokers.'
), file.path(out_dir, 'Missingness_notes.txt'))
cat('\nSaved to', normalizePath(out_dir), '\n')
