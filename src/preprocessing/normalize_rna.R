# =============================================================================
# Option B: library-size normalization of RNA only (TMM + log2-CPM)
# =============================================================================
# Keeps EVERYTHING else identical: same gene nodes, CpG nodes, edges, methylation,
# CNV and metadata. Only X_rna_gene_features.csv is replaced.
#
# Before: X_rna = log1p(raw STAR counts)            (no library-size correction)
# After : X_rna = log2-CPM with TMM factors          (edgeR::cpm, prior.count = 1)
#
# How it works (no Ensembl/Entrez re-query, so the gene set cannot change):
#   1. Raw counts of the 5,000+ exported genes are recovered exactly from the old
#      file: counts = expm1(log1p(counts)).
#   2. Library sizes and TMM factors come from ALL protein-coding genes of the same
#      samples in Session.RData (the RNA object Step_01 used), not only exported genes.
#   3. log2-CPM is computed per sample; no exposure or cross-patient outcome is used.
#
# Checks written to OUTPUT_DIR/rna_normalization_checks/:
#   - library sizes by tissue and tumor/normal depth ratio per patient,
#     and whether that ratio is associated with recorded exposure;
#   - agreement of old vs new values and of old vs new tumor-normal differences;
#   - overlap of the top-5,000 variance gene set (protein-coding Ensembl level)
#     selected on log1p raw counts vs on log2-CPM.
#
# Run (R, from Terminal):
#   Rscript normalize_rna_option_b.R SESSION_RDATA OLD_INPUT_DIR NEW_INPUT_DIR
# Example:
#   Rscript normalize_rna_option_b.R \
#     "$HOME/Library/CloudStorage/OneDrive-HarvardUniversity/GNN/TCGA/Lung/Session.RData" \
#     "$V/GNN_Input_Data" "$V/GNN_Input_Data_v2"
# Packages: BiocManager::install(c('edgeR', 'SummarizedExperiment'))

suppressPackageStartupMessages({
  library(edgeR)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('Usage: Rscript normalize_rna_option_b.R SESSION_RDATA OLD_INPUT_DIR NEW_INPUT_DIR')
session_file <- args[1]; old_dir <- args[2]; new_dir <- args[3]
if (!file.exists(session_file)) stop('Session file not found: ', session_file)
if (!dir.exists(old_dir)) stop('Old input folder not found: ', old_dir)
if (dir.exists(new_dir) && length(list.files(new_dir))) stop('NEW_INPUT_DIR must be new or empty: ', new_dir)
dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
check_dir <- file.path(new_dir, 'rna_normalization_checks')
dir.create(check_dir, showWarnings = FALSE)
TOP_VAR_GENES <- 5000   # same as Step_01

# =============================================================================
# 1. Old exported RNA (samples x genes, log1p raw counts)
# =============================================================================
old_rna <- read.csv(file.path(old_dir, 'X_rna_gene_features.csv'), row.names = 1, check.names = FALSE)
samples <- rownames(old_rna)
counts_selected <- t(expm1(as.matrix(old_rna)))                 # genes x samples, recovered counts
if (any(counts_selected < -1e-6)) stop('Negative recovered counts: is the old RNA file log1p of counts?')
counts_selected[counts_selected < 0] <- 0
cat('Exported genes:', nrow(counts_selected), ' samples:', ncol(counts_selected), '\n')

# =============================================================================
# 2. Library sizes and TMM factors from all protein-coding genes (Session.RData)
# =============================================================================
env <- new.env()
load(session_file, envir = env)
if (!exists('data_rna', envir = env)) stop('data_rna not found in ', session_file)
data_rna <- get('data_rna', envir = env); rm(env); invisible(gc())
missing <- setdiff(samples, colnames(data_rna))
if (length(missing)) stop(length(missing), ' exported samples are absent from data_rna, e.g. ', missing[1])
full_counts <- assay(data_rna)[, samples, drop = FALSE]        # first assay, as in Step_01
gene_type_col <- intersect(c('gene_type', 'gene_biotype'), colnames(rowData(data_rna)))
if (length(gene_type_col)) {
  pc <- rowData(data_rna)[[gene_type_col[1]]] == 'protein_coding'
  full_counts <- full_counts[which(pc), , drop = FALSE]
  cat('Protein-coding genes used for library size/TMM:', nrow(full_counts), '\n')
} else {
  cat('No gene_type column found: library size/TMM use all genes.\n')
}
full <- DGEList(counts = full_counts)
full <- calcNormFactors(full, method = 'TMM')
lib <- full$samples$lib.size; nf <- full$samples$norm.factors
names(lib) <- names(nf) <- samples

# =============================================================================
# 3. log2-CPM for the exported genes, using full-library sizes and TMM factors
# =============================================================================
sel <- DGEList(counts = counts_selected, lib.size = lib, norm.factors = nf)
logcpm <- cpm(sel, log = TRUE, prior.count = 1)                # genes x samples
new_rna <- t(logcpm)[samples, colnames(old_rna), drop = FALSE]
stopifnot(identical(dim(new_rna), dim(as.matrix(old_rna))), !anyNA(new_rna), all(is.finite(new_rna)))
write.csv(new_rna, file.path(new_dir, 'X_rna_gene_features.csv'), row.names = TRUE)

# Copy every other input file unchanged.
other <- setdiff(list.files(old_dir, full.names = TRUE, recursive = FALSE),
                 file.path(old_dir, 'X_rna_gene_features.csv'))
other <- other[!dir.exists(other)]
ok <- file.copy(other, new_dir, copy.date = TRUE)
if (!all(ok)) stop('Failed to copy: ', paste(basename(other[!ok]), collapse = ', '))
cat('Copied unchanged:', paste(basename(other), collapse = ', '), '\n')

# =============================================================================
# 4. Checks
# =============================================================================
meta <- read.csv(file.path(old_dir, 'metadata_samples.csv'), check.names = FALSE)
sid <- intersect(c('sample_id', 'sample_barcode'), names(meta))[1]
meta <- meta[match(samples, meta[[sid]]), ]
lib_table <- data.frame(sample_id = samples, patient_id = meta$patient_id, tissue_type = meta$tissue_type,
                        library_size = lib, tmm_factor = nf, effective_library = lib * nf,
                        smoking_intensity = meta$smoking_intensity)
write.csv(lib_table, file.path(check_dir, 'library_sizes.csv'), row.names = FALSE)

tp <- lib_table[lib_table$tissue_type == 'TP', ]; nt <- lib_table[lib_table$tissue_type == 'NT', ]
pairs <- merge(tp, nt, by = 'patient_id', suffixes = c('_tumor', '_normal'))
pairs$log2_depth_ratio <- log2(pairs$effective_library_tumor / pairs$effective_library_normal)
write.csv(pairs[, c('patient_id', 'effective_library_tumor', 'effective_library_normal', 'log2_depth_ratio', 'smoking_intensity_tumor')],
          file.path(check_dir, 'tumor_normal_depth_ratio.csv'), row.names = FALSE)

old_m <- as.matrix(old_rna)
gene_cor <- sapply(seq_len(ncol(old_m)), function(j) suppressWarnings(cor(old_m[, j], new_rna[, j], method = 'spearman')))
delta <- function(m) {
  m[match(pairs$sample_id_tumor, rownames(m)), , drop = FALSE] - m[match(pairs$sample_id_normal, rownames(m)), , drop = FALSE]
}
d_old <- delta(old_m); d_new <- delta(new_rna)
delta_cor <- sapply(seq_len(ncol(d_old)), function(j) suppressWarnings(cor(d_old[, j], d_new[, j], method = 'spearman')))

# Top-variance overlap at protein-coding Ensembl level (as close to Step_01 as possible
# without re-querying Ensembl/Entrez mappings).
full_logcpm <- cpm(full, log = TRUE, prior.count = 1)
var_old <- apply(log1p(full_counts), 1, var); var_new <- apply(full_logcpm, 1, var)
top_old <- names(sort(var_old, decreasing = TRUE))[seq_len(min(TOP_VAR_GENES, length(var_old)))]
top_new <- names(sort(var_new, decreasing = TRUE))[seq_len(min(TOP_VAR_GENES, length(var_new)))]
overlap <- length(intersect(top_old, top_new))

depth_rho <- suppressWarnings(cor(pairs$log2_depth_ratio, pairs$smoking_intensity_tumor, method = 'spearman'))
summary <- data.frame(
  check = c('Samples', 'Exported genes (unchanged)',
            'Median library size, tumor', 'Median library size, normal',
            'Median |log2 tumor/normal depth ratio| per patient',
            'Spearman: depth ratio vs recorded exposure (all patients in old input)',
            'Median per-gene Spearman, old vs new values',
            'Median per-gene Spearman, old vs new tumor-normal differences',
            paste0('Top-', TOP_VAR_GENES, ' variance genes shared, log1p-raw vs log2-CPM')),
  value = c(length(samples), ncol(old_m),
            format(median(tp$library_size), big.mark = ','), format(median(nt$library_size), big.mark = ','),
            round(median(abs(pairs$log2_depth_ratio)), 3), round(depth_rho, 3),
            round(median(gene_cor, na.rm = TRUE), 3), round(median(delta_cor, na.rm = TRUE), 3),
            paste0(overlap, ' / ', length(top_old), ' (', round(100 * overlap / length(top_old), 1), '%)')))
write.csv(summary, file.path(check_dir, 'normalization_summary.csv'), row.names = FALSE)
writeLines(c(
  'RNA library-size normalization (Option B)',
  'X_rna changed from log1p(raw STAR counts) to TMM-normalized log2-CPM (edgeR::cpm, prior.count = 1).',
  'Library sizes and TMM factors use all protein-coding genes of the same samples.',
  'Gene nodes, CpG nodes, edges, methylation, CNV and metadata files are copied unchanged.',
  'Normalization is per sample (TMM uses a cohort reference sample but no exposure or outcome information).',
  'Top-variance overlap is computed at protein-coding Ensembl level and approximates the Step_01 Entrez-level selection.'
), file.path(check_dir, 'README.txt'))
print(summary, row.names = FALSE)
cat('\nNew input folder:', normalizePath(new_dir), '\n')
