suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(SummarizedExperiment)
  library(tidyr)
  library(limma)
  library(GenomicRanges)
  library(matrixStats)
  library(biomaRt)
  library(impute)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(rtracklayer)
  library(AnnotationHub)
  library(STRINGdb)
  library(GenomeInfoDb)
})

# ==============================================================================
# USER SETTINGS
# ==============================================================================
# Usage: Rscript src/preprocessing/build_tcga_inputs.R WORK_DIR [OUTPUT_DIR]
#   WORK_DIR   folder containing Session.RData (data_rna, data_meth, data_cnv from GDC)
#   OUTPUT_DIR folder for the GNN inputs (default: WORK_DIR/GNN_Input_Data)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript build_tcga_inputs.R WORK_DIR [OUTPUT_DIR]")
setwd(args[1])

DATA_DIR          <- if (length(args) >= 2) args[2] else "GNN_Input_Data"
TOP_VAR_GENES     <- 5000    # high-variance gene pool (agnostic to exposure)
N_CPG_NODES       <- 500     # top tumor-vs-normal differential CpGs
STRING_SCORE      <- 700     # STRING confidence threshold (700 = high)

# Well-known smoking sentinel probes — added as prior anchors regardless of
# differential methylation status. Prior knowledge, not data-driven.
SMOKING_ANCHOR_CPGS <- c(
  "cg05575921", "cg03636183", "cg21161138",
  "cg19816560", "cg05951221", "cg06125433", "cg25648203"
)

# Smoking anchor genes (Entrez IDs): AHRR, CYP1A1, CYP1B1 and 5348
SMOKING_ANCHOR_GENES <- c("57491", "1543", "1545", "5348")
# Note: Entrez 5348 is FXYD1 (MYO1G is 64005). Kept as used in the published analysis.

TARGET_TISSUES <- tolower(c(
  "lower lobe, lung", "lung, nos", "main bronchus",
  "middle lobe, lung", "upper lobe, lung"
))

if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR)

# ==============================================================================
# HELPERS
# ==============================================================================
safe_num <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y), default, y)
}

# ==============================================================================
# 0) DOWNLOAD  (uncomment to run fresh; skip if Session.RData already exists)
# ==============================================================================
# projects <- c("TCGA-LUAD", "TCGA-LUSC")
#
# query_rna <- GDCquery(projects, "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification",
#                       workflow.type = "STAR - Counts")
# GDCdownload(query_rna, method = "client")
# data_rna <- GDCprepare(query_rna)
#
# query_meth <- GDCquery(projects, "DNA Methylation",
#                        data.type = "Methylation Beta Value",
#                        platform   = "Illumina Human Methylation 450")
# GDCdownload(query_meth, method = "client")
# data_meth <- GDCprepare(query_meth)
#
# query_cnv <- GDCquery(projects, "Copy Number Variation",
#                       data.type     = "Gene Level Copy Number",
#                       workflow.type = "ASCAT3")
# GDCdownload(query_cnv, method = "client")
# data_cnv <- GDCprepare(query_cnv)
#
# save(data_rna, data_meth, data_cnv, file = "Session.RData")

load("Session.RData")
stopifnot(exists("data_rna"), exists("data_meth"), exists("data_cnv"))

# ==============================================================================
# 1) FILTER BY TISSUE OF ORIGIN
# ==============================================================================
keep_tissue <- function(se) {
  tolower(colData(se)$tissue_or_organ_of_origin) %in% TARGET_TISSUES
}

data_rna  <- data_rna[,  keep_tissue(data_rna)]
data_meth <- data_meth[, keep_tissue(data_meth)]
data_cnv  <- data_cnv[,  keep_tissue(data_cnv)]

cat("After tissue filter — RNA:", ncol(data_rna),
    " Meth:", ncol(data_meth), " CNV:", ncol(data_cnv), "\n")

# ==============================================================================
# 2) CLINICAL METADATA
# ==============================================================================
demo <- bind_rows(
  GDCquery_clinic("TCGA-LUAD", "clinical"),
  GDCquery_clinic("TCGA-LUSC", "clinical")
) %>%
  dplyr::select(
    submitter_id,
    sex_at_birth,
    vital_status,
    age_at_index,
    tobacco_smoking_status,
    pack_years_smoked
  )

# ==============================================================================
# 3) PAIRED SAMPLES — one tumor (TP) + one normal (NT) per patient
#    Uses RNA barcodes as the master sample list.
#    Methylation and CNV are aligned to this list in section 5.
# ==============================================================================
all_barcodes    <- colData(data_rna)$barcode
paired_barcodes <- TCGAquery_MatchedCoupledSampleTypes(
  barcode = all_barcodes, types = c("TP", "NT")
)

rna_paired <- data_rna[, paired_barcodes]
colData(rna_paired)$submitter_id <- substr(colData(rna_paired)$barcode, 1, 12)

# Attach clinical metadata
cd <- as.data.frame(colData(rna_paired)) %>%
  left_join(demo, by = "submitter_id")
rownames(cd) <- cd$barcode
colData(rna_paired) <- S4Vectors::DataFrame(cd)

# Keep exactly 1 TP + 1 NT per patient
final_barcodes <- as.data.frame(colData(rna_paired)) %>%
  transmute(
    barcode      = as.character(barcode),
    submitter_id = as.character(submitter_id),
    sample_code  = substr(as.character(barcode), 14, 15)
  ) %>%
  dplyr::filter(sample_code %in% c("01", "11")) %>%
  dplyr::group_by(submitter_id, sample_code) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pull(barcode)

rna_gnn     <- rna_paired[, final_barcodes]
rna_samples <- colnames(rna_gnn)
cat("Paired samples:", length(rna_samples), "\n")

# ==============================================================================
# 4) EXPOSURE + METADATA TABLE
# ==============================================================================
smoking_col   <- grep("smoking_status", colnames(colData(rna_gnn)), value = TRUE)[1]
pack_year_col <- grep("pack_years",     colnames(colData(rna_gnn)), value = TRUE)[1]
age_col       <- grep("age_at_index",   colnames(colData(rna_gnn)), value = TRUE)[1]
gender_col    <- grep("sex_at_birth",         colnames(colData(rna_gnn)), value = TRUE)[1]
vital_col     <- grep("vital_status",   colnames(colData(rna_gnn)), value = TRUE)[1]

if (is.na(smoking_col)) stop("tobacco_smoking_status not found in colData.")

metadata <- as.data.frame(colData(rna_gnn)) %>%
  mutate(sample_barcode = rownames(.)) %>%
  transmute(
    sample_barcode,
    patient_id        = substr(sample_barcode, 1, 12),
    tissue_type       = ifelse(as.numeric(substr(sample_barcode, 14, 15)) < 10,
                               "TP", "NT"),
    smoking_status    = .data[[smoking_col]],
    smoking_rank      = case_when(
      smoking_status == "Lifelong Non-Smoker"                        ~ 1,
      smoking_status == "Current Reformed Smoker for > 15 yrs"       ~ 2,
      smoking_status == "Current Reformed Smoker for < or = 15 yrs"  ~ 3,
      smoking_status == "Current Smoker"                             ~ 4,
      TRUE                                                           ~ 2.5
    ),
    # Primary exposure variable (pack_years); 0 for never-smokers
    smoking_intensity = coalesce(safe_num(.data[[pack_year_col]]), 0),
    age               = coalesce(safe_num(.data[[age_col]]),      65),
    gender            = .data[[gender_col]],
    vital_status      = .data[[vital_col]],
    diagnosis         = as.integer(tissue_type == "TP")
  )

cat("Metadata rows:", nrow(metadata), "\n")

# ==============================================================================
# 5) ALIGN METHYLATION + CNV TO RNA SAMPLE ORDER (16-character barcode key)
# ==============================================================================
rna_key16 <- substr(rna_samples, 1, 16)

# --- Methylation ---
meth_key16   <- substr(colnames(data_meth), 1, 16)
beta_all     <- assay(data_meth)
meth_idx     <- match(rna_key16, meth_key16)
beta_aligned <- beta_all[, meth_idx, drop = FALSE]
colnames(beta_aligned) <- rna_samples

missing_meth <- which(is.na(meth_idx))
if (length(missing_meth) > 0) {
  normal_cols      <- which(grepl("-11", meth_key16))
  meth_baseline    <- apply(beta_all[, normal_cols, drop = FALSE], 1,
                            median, na.rm = TRUE)
  beta_aligned[, missing_meth] <- meth_baseline
  cat("Meth: imputed", length(missing_meth),
      "missing samples with normal-tissue median.\n")
}

# --- CNV ---
# Raw values from ASCAT3 are integer absolute copy numbers (0, 1, 2, 3 ...).
# Normal samples are uniformly diploid (all 2), giving IQR = 0 and collapsing
# the signal when normal+tumor are combined.
# Fix: impute NAs as 2, center at diploid, then sign*log2-compress.
#   sign(x) * log2(|x| + 1)  maps:
#     diploid (0 after centering) -> 0.0   (exact zero)
#     +1 copy gain               -> +1.0
#     +2 copy gain               -> +1.58
#     +56 copy amplification     -> +5.84  (not 56x — compressed)
#     -1 copy deletion           -> -1.0
#     -2 copy homozygous del     -> -1.58
# This gives informative IQRs for CNA-frequent genes while keeping
# normal samples exactly zero.
cnv_key16   <- substr(colnames(data_cnv), 1, 16)
cnv_all     <- tryCatch(assay(data_cnv, "copy_number"), error = function(e) assay(data_cnv))
cnv_idx     <- match(rna_key16, cnv_key16)
cnv_raw     <- cnv_all[, cnv_idx, drop = FALSE]
colnames(cnv_raw) <- rna_samples

missing_cnv <- which(is.na(cnv_idx))
if (length(missing_cnv) > 0) {
  cnv_raw[, missing_cnv] <- 2L
  cat("CNV: imputed", length(missing_cnv), "missing samples as diploid (2).\n")
}

# Impute within-gene NAs as 2 (diploid neutral for unmapped regions)
na_count <- sum(is.na(cnv_raw))
if (na_count > 0) {
  cnv_raw[is.na(cnv_raw)] <- 2L
  cat(sprintf("CNV: imputed %d within-matrix NAs as diploid.\n", na_count))
}

# Strip Ensembl version suffixes from row names so they match final_mapping
# which uses base IDs (e.g. ENSG00000000003, not ENSG00000000003.15).
# Without this, the intersect() in section 10 returns zero rows and X_cnv
# is filled entirely with the diploid-neutral default (0), losing all signal.
rownames(cnv_raw) <- gsub("\\..*", "", rownames(cnv_raw))

# Center then sign-log2 compress
cnv_centered <- matrix(as.numeric(cnv_raw) - 2,
                       nrow = nrow(cnv_raw), ncol = ncol(cnv_raw),
                       dimnames = dimnames(cnv_raw))
cnv_aligned  <- sign(cnv_centered) * log2(abs(cnv_centered) + 1)

# Diagnostic summary
normal_samples <- colnames(cnv_aligned)[grepl("-11[A-Z]-", colnames(cnv_aligned))]
if (length(normal_samples) > 0) {
  nr <- range(cnv_aligned[, normal_samples])
  cat(sprintf("CNV normal sample range after transform: [%.2f, %.2f]%s\n",
              nr[1], nr[2],
              ifelse(all(cnv_aligned[, normal_samples] == 0),
                     " (all zero — correct)", " (WARNING: non-zero)")))
}
tumor_samples <- colnames(cnv_aligned)[grepl("-01[A-Z]-", colnames(cnv_aligned))]
if (length(tumor_samples) > 0) {
  iqr_t <- apply(cnv_aligned[, tumor_samples, drop = FALSE], 1, IQR)
  cat(sprintf("CNV tumor-only IQR — median=%.3f  max=%.3f  genes>0.1: %d/%d\n",
              median(iqr_t), max(iqr_t),
              sum(iqr_t > 0.1), length(iqr_t)))
}

# ==============================================================================
# 6) CpG NODE SELECTION — tumor-vs-normal differential methylation
#
#    Model: ~ tissue_type  (NO smoking term — exposure-agnostic selection)
#    Selection: top N_CPG_NODES by |logFC| among FDR < 0.05
#    Rationale: CpGs selected on cancer relevance only. The GNN then
#    discovers which cancer-relevant CpGs are additionally modulated by
#    smoking. Avoids circularity with the exposure variable.
#
#    Known smoking sentinel probes added as hard prior anchors.
# ==============================================================================
tissue_type <- factor(metadata$tissue_type[
  match(colnames(beta_aligned), metadata$sample_barcode)
], levels = c("NT", "TP"))

design_cpg <- model.matrix(~ tissue_type)
fit_cpg    <- eBayes(lmFit(beta_aligned, design_cpg))
tt_cpg     <- topTable(fit_cpg, coef = 2, number = Inf,
                       sort.by = "P", adjust.method = "fdr")

sig_cpgs <- tt_cpg %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  slice_max(abs(logFC), n = N_CPG_NODES)

key_cpg_ids <- unique(c(rownames(sig_cpgs), SMOKING_ANCHOR_CPGS))
key_cpg_ids <- intersect(key_cpg_ids, rownames(data_meth))

cat("CpG nodes selected (tumor-vs-normal, FDR<0.05):", length(key_cpg_ids), "\n")

# ==============================================================================
# 7) CpG → NEAREST GENE EDGES  (hg19 CpG positions lifted to hg38)
# ==============================================================================
common_ensembl      <- intersect(rownames(rna_gnn), rownames(cnv_all))
common_ensembl_base <- gsub("\\..*", "", common_ensembl)

# Gene TSS positions in hg38
ensembl38 <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl",
                                 mirror = "useast")
gene_pos <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "chromosome_name",
                 "start_position", "end_position", "strand"),
  filters    = "ensembl_gene_id",
  values     = common_ensembl_base,
  mart       = ensembl38
) %>%
  dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

gene_tss <- GRanges(
  seqnames = paste0("chr", gene_pos$chromosome_name),
  ranges   = IRanges(
    start = ifelse(gene_pos$strand == 1,
                   gene_pos$start_position,
                   gene_pos$end_position),
    width = 1
  ),
  strand = ifelse(gene_pos$strand == 1, "+", "-")
)
names(gene_tss)           <- gene_pos$ensembl_gene_id
GenomeInfoDb::genome(gene_tss) <- "hg38"
cat("Gene TSS built (hg38):", length(gene_tss), "\n")

# LiftOver CpG positions hg19 -> hg38
cpg_ranges_hg19 <- rowRanges(data_meth)[key_cpg_ids]
seqlevelsStyle(cpg_ranges_hg19) <- "UCSC"
GenomeInfoDb::genome(cpg_ranges_hg19) <- "hg19"

ah      <- AnnotationHub()
chain_q <- query(ah, c("hg19ToHg38", "chain"))
if (length(chain_q) == 0) stop("hg19ToHg38 chain not found in AnnotationHub.")
chain <- chain_q[[1]]

cpg_lift <- liftOver(cpg_ranges_hg19, chain)
one_map  <- elementNROWS(cpg_lift) == 1
cpg_hg38 <- unlist(cpg_lift[one_map])
names(cpg_hg38) <- names(cpg_lift)[one_map]
GenomeInfoDb::genome(cpg_hg38) <- "hg38"

dropped <- setdiff(key_cpg_ids, names(cpg_hg38))
if (length(dropped) > 0)
  cat("LiftOver dropped", length(dropped), "CpGs (multi-map or unmapped).\n")

key_cpg_ids <- names(cpg_hg38)
cat("CpGs after liftOver:", length(key_cpg_ids), "\n")

# Nearest gene per CpG
nearest_hits <- distanceToNearest(cpg_hg38, gene_tss)

edge_cpg_gene_raw <- data.frame(
  cpg_id  = names(cpg_hg38)[queryHits(nearest_hits)],
  gene_id = names(gene_tss)[subjectHits(nearest_hits)],
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct(cpg_id, gene_id)

cat("CpG->Gene edges (raw, hg38):", nrow(edge_cpg_gene_raw), "\n")

# ==============================================================================
# 8) ENSEMBL -> ENTREZ MAPPING  (protein-coding genes only)
# ==============================================================================
ensembl_mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

rna_mat     <- assay(rna_gnn)[common_ensembl, , drop = FALSE]
raw_ensembl <- gsub("\\..*", "", rownames(rna_mat))

mapping_table <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id",
                 "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = raw_ensembl,
  mart       = ensembl_mart
) %>%
  dplyr::filter(!is.na(entrezgene_id), gene_biotype == "protein_coding") %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(entrez_id = as.character(entrezgene_id))

# Collapse RNA matrix from Ensembl to Entrez IDs (mean over duplicates)
mat_rna_entrez <- as.data.frame(rna_mat) %>%
  mutate(ensembl_id = raw_ensembl) %>%
  inner_join(mapping_table, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::select(-ensembl_id, -entrezgene_id, -external_gene_name, -gene_biotype) %>%
  group_by(entrez_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("entrez_id") %>%
  as.matrix()

cat("RNA matrix (Entrez):", paste(dim(mat_rna_entrez), collapse = " x "), "\n")

# ==============================================================================
# 9) GENE NODE SELECTION
#
#    Three sources, all agnostic to the exposure variable:
#    (a) Top TOP_VAR_GENES by expression variance across all samples
#    (b) Target genes of selected CpG nodes (genes the CpGs point to)
#    (c) Known smoking-anchor genes (prior knowledge)
#
#    Final = union(a, b, c) restricted to protein-coding genes in RNA matrix.
# ==============================================================================

# (a) High-variance genes (variance computed on log1p scale)
gene_var      <- rowVars(log1p(mat_rna_entrez))
high_var_pool <- head(names(sort(gene_var, decreasing = TRUE)), TOP_VAR_GENES)

# (b) CpG target genes — map Ensembl -> Entrez
edge_cpg_gene_raw$gene_id_base <- gsub("\\..*", "", edge_cpg_gene_raw$gene_id)

edge_cpg_entrez <- edge_cpg_gene_raw %>%
  inner_join(mapping_table, by = c("gene_id_base" = "ensembl_gene_id")) %>%
  transmute(cpg_id, entrez_id = as.character(entrez_id)) %>%
  distinct()

genes_from_cpg <- unique(edge_cpg_entrez$entrez_id)

# (c) Union
final_gene_set <- unique(c(SMOKING_ANCHOR_GENES, genes_from_cpg, high_var_pool))
final_gene_set <- intersect(final_gene_set, rownames(mat_rna_entrez))

cat("Final gene nodes:", length(final_gene_set), "\n")
cat("  from high-var pool  :", sum(final_gene_set %in% high_var_pool), "\n")
cat("  from CpG targets    :", sum(final_gene_set %in% genes_from_cpg), "\n")
cat("  smoking anchors     :", sum(final_gene_set %in% SMOKING_ANCHOR_GENES), "\n")

final_mapping <- mapping_table %>%
  dplyr::filter(entrez_id %in% final_gene_set) %>%
  distinct(entrez_id, .keep_all = TRUE)

# ==============================================================================
# 10) GENE NODE FEATURE MATRICES
#
#     Each gene node carries 3 features: [RNA, mean_promoter_beta, CNV].
#     All three are patient-specific.
#
#     RNA   : log1p of raw counts — compresses dynamic range,
#             makes IQRs comparable to meth and cnv units.
#     Meth  : mean beta across TSS1500/TSS200/5'UTR/1stExon probes per gene.
#             Beta is already in [0,1]; Step_02 applies robust scaling.
#     CNV   : sign*log2(|centered|+1) from section 5.
#             0 = diploid, positive = gain, negative = deletion.
#             Step_02 applies robust scaling with cnv_min_scale=0.05.
# ==============================================================================

# --- RNA (log1p) ---
X_rna <- log1p(mat_rna_entrez[final_gene_set, , drop = FALSE])
# ============================================================
# Mean promoter methylation per gene
# FIX: use UCSC_RefGene_Name gene SYMBOLS, not Accession.
# Then map gene symbol -> Entrez ID using mapping_table.
# ============================================================

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_df  <- as.data.frame(ann450k)

promoter_groups <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

promoter_map <- ann_df %>%
  dplyr::select(
    probe_id      = Name,
    refgene_name  = UCSC_RefGene_Name,
    refgene_group = UCSC_RefGene_Group
  ) %>%
  dplyr::mutate(
    probe_id      = as.character(probe_id),
    refgene_name  = as.character(refgene_name),
    refgene_group = as.character(refgene_group)
  ) %>%
  tidyr::separate_rows(refgene_name, refgene_group, sep = ";") %>%
  dplyr::filter(
    !is.na(probe_id), probe_id != "",
    !is.na(refgene_name), refgene_name != "",
    !is.na(refgene_group), refgene_group != "",
    refgene_group %in% promoter_groups
  ) %>%
  dplyr::distinct()

promoter_to_entrez <- promoter_map %>%
  dplyr::inner_join(
    mapping_table %>%
      dplyr::select(external_gene_name, entrez_id) %>%
      dplyr::mutate(
        external_gene_name = as.character(external_gene_name),
        entrez_id          = as.character(entrez_id)
      ) %>%
      dplyr::filter(
        !is.na(external_gene_name), external_gene_name != "",
        !is.na(entrez_id), entrez_id != ""
      ) %>%
      dplyr::distinct(),
    by = c("refgene_name" = "external_gene_name")
  ) %>%
  dplyr::transmute(
    probe_id  = as.character(probe_id),
    entrez_id = as.character(entrez_id)
  ) %>%
  dplyr::filter(
    probe_id %in% rownames(beta_aligned),
    entrez_id %in% final_gene_set
  ) %>%
  dplyr::distinct()

cat(sprintf("Promoter probes mapped: %d\n", nrow(promoter_to_entrez)))
cat(sprintf("Promoter genes mapped : %d\n",
            length(unique(promoter_to_entrez$entrez_id))))

# Aggregate promoter beta values: probes -> genes
meth_promoter <- as.data.frame(
  beta_aligned[promoter_to_entrez$probe_id, , drop = FALSE]
) %>%
  tibble::rownames_to_column("probe_id") %>%
  dplyr::inner_join(promoter_to_entrez, by = "probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::group_by(entrez_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames("entrez_id")

# Align to final gene set and RNA sample order
X_meth_gene <- matrix(
  NA_real_,
  nrow = length(final_gene_set),
  ncol = length(rna_samples),
  dimnames = list(as.character(final_gene_set), as.character(rna_samples))
)

shared_samples_meth <- intersect(colnames(meth_promoter), rna_samples)
common_meth <- intersect(rownames(meth_promoter), rownames(X_meth_gene))

X_meth_gene[common_meth, shared_samples_meth] <- as.matrix(
  meth_promoter[common_meth, shared_samples_meth, drop = FALSE]
)

# Per-sample median imputation for missing gene methylation values
sample_medians <- apply(X_meth_gene, 2, median, na.rm = TRUE)

for (j in seq_len(ncol(X_meth_gene))) {
  miss <- is.na(X_meth_gene[, j])
  if (any(miss)) {
    X_meth_gene[miss, j] <- sample_medians[j]
  }
}

# Safety fallback if any column median was NA
X_meth_gene[is.na(X_meth_gene)] <- median(X_meth_gene, na.rm = TRUE)

# Sanity checks
n_covered <- length(common_meth)
n_total   <- nrow(X_meth_gene)

cat(sprintf("Genes with promoter methylation coverage: %d / %d (%.1f%%)\n",
            n_covered, n_total, 100 * n_covered / n_total))

cat(sprintf("Range after imputation: [%.4f, %.4f]\n",
            min(X_meth_gene), max(X_meth_gene)))

cat(sprintf("Median after imputation: %.4f\n",
            median(X_meth_gene)))

gene_sd <- apply(X_meth_gene, 1, sd)
cat(sprintf("Median per-gene SD: %.6f\n", median(gene_sd)))
cat(sprintf("IQR per-gene SD   : %.6f\n", IQR(gene_sd)))
cat(sprintf("Genes with SD > 0.001: %d / %d\n",
            sum(gene_sd > 0.001), length(gene_sd)))

if (n_covered < 0.10 * n_total) {
  warning(sprintf(
    "Only %.1f%% of gene nodes have promoter methylation coverage. Check probe-to-gene mapping.",
    100 * n_covered / n_total
  ))
}

# Sanity checks
stopifnot("NAs remain after meth imputation" = sum(is.na(X_meth_gene)) == 0)
cat(sprintf("  Range after imputation       : [%.4f, %.4f]\n",
            min(X_meth_gene), max(X_meth_gene)))
cat(sprintf("  Median after imputation      : %.4f\n",
            median(X_meth_gene)))
cat(sprintf("  IQR of per-gene SD           : %.6f  (should be > 0)\n",
            IQR(apply(X_meth_gene, 1, sd))))

# --- CNV (sign-log2, already computed in section 5) ---
# Collapse from Ensembl to Entrez (mean over duplicates), then align
relevant_ensembl <- final_mapping$ensembl_gene_id
mat_cnv <- cnv_aligned[intersect(rownames(cnv_aligned), relevant_ensembl), , drop = FALSE]

X_cnv_collapsed <- as.data.frame(mat_cnv) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  mutate(ensembl_id = gsub("\\..*", "", ensembl_id)) %>%
  inner_join(final_mapping, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::select(entrez_id, where(is.numeric)) %>%
  group_by(entrez_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("entrez_id")

# Align to master gene list; genes without CNV data are set to 0 (diploid neutral)
X_cnv <- matrix(
  0,
  nrow     = length(final_gene_set),
  ncol     = length(rna_samples),
  dimnames = list(final_gene_set, rna_samples)
)
# Reorder X_cnv_collapsed columns to match rna_samples exactly.
# After dplyr summarise the column order may differ, causing the
# 'not a multiple of replacement length' error on assignment.
shared_samples_cnv <- intersect(colnames(X_cnv_collapsed), rna_samples)
X_cnv_collapsed    <- X_cnv_collapsed[, shared_samples_cnv, drop = FALSE]

common_cnv <- intersect(rownames(X_cnv_collapsed), rownames(X_cnv))
X_cnv[common_cnv, shared_samples_cnv] <- as.matrix(
  X_cnv_collapsed[common_cnv, , drop = FALSE]
)

cat("CNV matrix aligned:", nrow(X_cnv), "genes\n")

# Sanity check: confirm X_rna, X_meth_gene, X_cnv have identical dimensions
stopifnot(
  identical(dim(X_rna), dim(X_meth_gene)),
  identical(dim(X_rna), dim(X_cnv)),
  identical(colnames(X_rna), rna_samples)
)

# ==============================================================================
# 11) CpG NODE FEATURE MATRIX
#     One feature per CpG node: raw beta value (patient-specific).
#     Only CpGs that successfully lifted to hg38 and have valid edges.
# ==============================================================================
sentinel_probes <- intersect(unique(edge_cpg_entrez$cpg_id), rownames(beta_aligned))
X_cpg           <- beta_aligned[sentinel_probes, , drop = FALSE]

if (any(is.na(X_cpg))) {
  X_cpg <- impute.knn(as.matrix(X_cpg), k = 10)$data
  cat("CpG: KNN-imputed missing beta values.\n")
}

cat("CpG nodes:", nrow(X_cpg), "\n")

# ==============================================================================
# 12) NODE INDEX MAPS
#     Gene nodes: index 0 .. (n_genes - 1)        node_type = 1
#     CpG  nodes: index n_genes .. (n_total - 1)   node_type = 0
#     These indices are used directly by the PyG edge_index tensors.
# ==============================================================================
gene_nodes_map <- data.frame(
  entrez_id = rownames(X_rna),
  node_idx  = seq_len(nrow(X_rna)) - 1L,
  node_type = 1L
)

cpg_nodes_map <- data.frame(
  cpg_id    = rownames(X_cpg),
  node_idx  = nrow(X_rna) + seq_len(nrow(X_cpg)) - 1L,
  node_type = 0L
)

cat("Node counts — genes:", nrow(gene_nodes_map),
    " CpGs:", nrow(cpg_nodes_map),
    " total:", nrow(gene_nodes_map) + nrow(cpg_nodes_map), "\n")

# ==============================================================================
# 13) EDGES
#
#     PPI  : gene <-> gene  (undirected, STRING high-confidence)
#             Both directions stored explicitly.
#     CpG  : CpG  ->  gene  (directed, nearest-gene assignment)
#             NOTE: Step_02 automatically adds the reverse Gene->CpG direction
#             so CpG nodes receive message-passing context from their gene.
#             Only the CpG->Gene direction is stored here.
#
#     All edges: weight = 1.0  (binary structural prior).
#     Patient-specific rewiring lives in GAT attention scores, not weights.
# ==============================================================================
# library(httr)
# # Set global connection parameters to be more patient
# set_config(config(connecttimeout = 60, timeout = 300))

options(timeout = 600) 

# --- PPI edges ---
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = STRING_SCORE)

mapped_genes <- string_db$map(
  data.frame(entrez_id = rownames(X_rna)),
  "entrez_id", removeUnmappedRows = TRUE
)

id_bridge <- mapped_genes %>%
  dplyr::select(STRING_id, entrez_id) %>%
  group_by(STRING_id) %>%
  dplyr::slice(1) %>%
  ungroup()

cat("STRING mapped:", nrow(id_bridge), "/", nrow(X_rna), "genes\n")

interactions <- string_db$get_interactions(id_bridge$STRING_id)

ppi_edges_raw <- interactions %>%
  inner_join(id_bridge, by = c("from" = "STRING_id")) %>%
  dplyr::rename(from_entrez = entrez_id) %>%
  inner_join(id_bridge, by = c("to"   = "STRING_id")) %>%
  dplyr::rename(to_entrez = entrez_id) %>%
  inner_join(gene_nodes_map, by = c("from_entrez" = "entrez_id")) %>%
  dplyr::rename(from_idx = node_idx) %>%
  inner_join(gene_nodes_map, by = c("to_entrez"   = "entrez_id")) %>%
  dplyr::rename(to_idx = node_idx) %>%
  dplyr::select(from_idx, to_idx) %>%
  distinct()

# Undirected: store both directions
edges_ppi <- bind_rows(
  ppi_edges_raw,
  ppi_edges_raw %>% transmute(from_idx = to_idx, to_idx = from_idx)
) %>%
  mutate(weight = 1.0) %>%
  distinct()

cat("PPI edges (undirected):", nrow(edges_ppi), "\n")

# --- CpG -> Gene edges ---
edges_cpg_gene <- edge_cpg_entrez %>%
  dplyr::filter(
    cpg_id    %in% cpg_nodes_map$cpg_id,
    entrez_id %in% gene_nodes_map$entrez_id
  ) %>%
  inner_join(cpg_nodes_map,  by = "cpg_id") %>%
  dplyr::rename(from_idx = node_idx) %>%
  inner_join(gene_nodes_map, by = "entrez_id") %>%
  dplyr::rename(to_idx = node_idx) %>%
  transmute(from_idx, to_idx, weight = 1.0, cpg_id, entrez_id) %>%
  distinct()

cat("CpG->Gene edges:", nrow(edges_cpg_gene), "\n")

# ==============================================================================
# 14) METADATA — aligned to master sample order (colnames of X_rna)
# ==============================================================================
final_metadata <- metadata %>%
  dplyr::filter(sample_barcode %in% colnames(X_rna))

final_metadata <- final_metadata[
  match(colnames(X_rna), final_metadata$sample_barcode), ]

stopifnot(identical(final_metadata$sample_barcode, colnames(X_rna)))

map_sample_index <- data.frame(
  sample_id    = colnames(X_rna),
  sample_index = seq_along(colnames(X_rna)) - 1L
)

final_metadata <- left_join(
  map_sample_index, final_metadata,
  by = c("sample_id" = "sample_barcode")
)

# Rename sample_id -> sample_id for consistency with Step_02 default
# (sample_id_col = "sample_id" in LungMultiOmicsDataset)
names(final_metadata)[names(final_metadata) == "sample_id"] <- "sample_id"

# ==============================================================================
# 15) EXPORT
#
#     All feature matrices written as SAMPLES x FEATURES (rows = samples).
#     Python (Step_02) reads them with index_col=0 giving samples as index.
#
#     X_rna_gene_features.csv   [S x G]  log1p RNA expression
#     X_meth_gene_features.csv  [S x G]  mean promoter beta per gene [0,1]
#     X_cnv_gene_features.csv   [S x G]  sign*log2(|copy_number - 2| + 1)
#     X_cpg_features.csv        [S x C]  raw beta per CpG node [0,1]
#     edges_ppi.csv             from_idx, to_idx, weight
#     edges_cpg_to_gene.csv     from_idx (CpG), to_idx (gene), weight, ids
#     map_gene_nodes.csv        entrez_id, node_idx, node_type=1
#     map_cpg_nodes.csv         cpg_id,    node_idx, node_type=0
#     metadata_samples.csv      clinical + exposure variables
#     map_sample_index.csv      sample barcode -> integer index
# ==============================================================================
write.csv(t(X_rna),         file.path(DATA_DIR, "X_rna_gene_features.csv"),  row.names = TRUE)
write.csv(t(X_meth_gene),   file.path(DATA_DIR, "X_meth_gene_features.csv"), row.names = TRUE)
write.csv(t(X_cnv),         file.path(DATA_DIR, "X_cnv_gene_features.csv"),  row.names = TRUE)
write.csv(t(X_cpg),         file.path(DATA_DIR, "X_cpg_features.csv"),       row.names = TRUE)

write.csv(edges_ppi,        file.path(DATA_DIR, "edges_ppi.csv"),            row.names = FALSE)
write.csv(edges_cpg_gene,   file.path(DATA_DIR, "edges_cpg_to_gene.csv"),    row.names = FALSE)

write.csv(gene_nodes_map,   file.path(DATA_DIR, "map_gene_nodes.csv"),       row.names = FALSE)
write.csv(cpg_nodes_map,    file.path(DATA_DIR, "map_cpg_nodes.csv"),        row.names = FALSE)
write.csv(final_metadata,   file.path(DATA_DIR, "metadata_samples.csv"),     row.names = FALSE)
write.csv(map_sample_index, file.path(DATA_DIR, "map_sample_index.csv"),     row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Genes    :", nrow(X_rna), "\n")
cat("CpGs     :", nrow(X_cpg), "\n")
cat("Samples  :", ncol(X_rna), "\n")
cat("PPI edges:", nrow(edges_ppi), "\n")
cat("CpG->Gene:", nrow(edges_cpg_gene), "\n")

