# =============================================================================
# ExposureGAT simulation: Supplementary Figures 3-5 (node scores by class, paired-transition
# magnitude by class, exposure association of edge attention changes). Run per manifest of datasets.
# =============================================================================
# Run from Terminal:
# Rscript plot_all_results.R DATA_DIR ANALYSIS_DIR OUTPUT_DIR [MANIFEST.csv]
# Optional manifest columns: seed,data_dir,analysis_dir (absolute paths).
# Use outputs from the edge_transition_extension version of analyze.py.
# No fitting, p-values, or changes to model scores occur in this script.
# install.packages(c('dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite'))

#Rscript plot_all_results.R \
#"results/seed_1/data" \
#"results/seed_1/analysis" \
#"results/R_plots" \
#"results/manifest.csv"

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(ggplot2)
library(ggpubr)
library(dplyr)
# =============================================================================
# Settings: edit these paths when running interactively in RStudio
# =============================================================================
data_dir <- NULL; analysis_dir <- NULL; output_dir <- NULL; manifest_file <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  data_dir <- args[1]
  analysis_dir <- args[2]
  output_dir <- args[3]
}
if (length(args) >= 4) manifest_file <- args[4]
if (length(args) < 4) stop('Usage: Rscript plot_simulation_supplement.R DATA_DIR ANALYSIS_DIR OUTPUT_DIR MANIFEST.csv')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
# Figure 2 : Transition score vs smoking
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
save_plot(p_transition, 'Patient_transition_vs_normalized_exposure', width = 5, height = 4)

# =============================================================================
# Figure 2: Predicted smoking vs actual smoking: matched rank units
# =============================================================================
# predicted_rho <- rho(example_patients$predicted_exposure_rank_z, example_patients$actual_exposure_rank_z)
# p_predicted <- ggplot(example_patients, aes(actual_exposure_rank_z, predicted_exposure_rank_z)) +
#   geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
#   geom_point(alpha = 0.7, color = '#8361BB') +
#   geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black') +
#   labs(title = 'Predicted exposure rank vs actual exposure rank',
#        subtitle = paste0('Dataset ', example_seed, '; Spearman rho = ', round(predicted_rho, 3)),
#        x = 'Actual standardized exposure rank', y = 'Predicted standardized exposure-rank score',
#        caption = 'In-sample diagnostic; dashed line = agreement. Predictions are not pack-years.')
# rank_limits <- range(c(example_patients$actual_exposure_rank_z, example_patients$predicted_exposure_rank_z))
# p_predicted <- p_predicted + coord_equal(xlim = rank_limits, ylim = rank_limits)
# save_plot(p_predicted, 'Predicted_vs_actual_smoking_rank')
# 
# p_raw <- ggplot(example_patients, aes(actual_exposure, predicted_exposure_rank_z)) +
#   geom_point(alpha = 0.7, color = '#8361BB') +
#   labs(title = 'Fitted exposure score vs recorded smoking',
#        subtitle = paste0('Spearman rho = ', round(predicted_rho, 3)),
#        x = paste0('Recorded exposure: ', exposure_column), y = 'Predicted standardized exposure-rank score',
#        caption = 'Different units on the two axes; no raw smoking-amount prediction is made.')
# save_plot(p_raw, 'Predicted_rank_vs_recorded_smoking')
# 
# 
# 
# ################################################################################
# #### Figure 2: NODE-level ROC comparisons across independent simulation datasets.
# ################################################################################
#manifest_file <- 'results/manifest.csv'
plot_dir <- output_dir
manifest <- read_csv(manifest_file,show_col_types=FALSE)
node_rows = list()
for(i in seq_len(nrow(manifest))) {
  truth <- read_csv(file.path(manifest$data_dir[i],'sim_ground_truth_nodes.csv'),show_col_types=FALSE) %>%
    select(node_idx,set_membership)
  scores <- read_csv(file.path(manifest$analysis_dir[i],'cohort_summary.csv'),show_col_types=FALSE) %>%
    select(node_idx,median_patient_rank)
  if(anyDuplicated(truth$node_idx)||anyDuplicated(scores$node_idx)) stop('Duplicate node indices.')
  dat <- scores %>% left_join(truth,by='node_idx')
  if(anyNA(dat$set_membership)||any(!dat$set_membership %in% c('A','B','C'))) stop('Incomplete A/B/C truth mapping.')
  if(any(!is.finite(dat$median_patient_rank))) stop('Nonfinite node ranks.')
  if(nrow(dat)<2) stop('Not enough nodes.')
  node_rows[[i]] <- dat %>% mutate(NodeScore=(n()-median_patient_rank)/(n()-1),seed=manifest$seed[i])
}
all_nodes <- bind_rows(node_rows)
write_csv(all_nodes,file.path(plot_dir,'Table_Node_Scores_Used_For_ROC.csv'))
figure_theme <- theme_bw(base_size=14)+theme(plot.title=element_text(face='bold'))
save_plot <- function(p,name,width=8,height=5) {
  ggsave(file.path(plot_dir,paste0(name,'.pdf')),p,width=width,height=height)
  ggsave(file.path(plot_dir,paste0(name,'.png')),p,width=width,height=height,dpi=300)
}
# Tied thresholds grouped; fixed direction; exact AUROC per dataset.
roc_comparisons <- list('A vs B'=c('A','B'),'C vs B'=c('C','B'),'A vs C'=c('A','C'))
roc_rows <- list(); auc_rows <- list(); roc_grid_rows <- list()
fpr_grid <- seq(0,1,length.out=501)
for(seed_id in manifest$seed) {
  for(comparison in names(roc_comparisons)) {
    classes <- roc_comparisons[[comparison]]
    dat <- all_nodes %>% filter(seed==seed_id,set_membership %in% classes,is.finite(NodeScore)) %>%
      mutate(positive=set_membership==classes[1])
    n_pos <- sum(dat$positive); n_neg <- sum(!dat$positive)
    if(n_pos==0 || n_neg==0) stop('ROC requires both classes: ',comparison,' seed ',seed_id)
    grouped <- dat %>% group_by(NodeScore) %>%
      summarise(tp=sum(positive),fp=sum(!positive),.groups='drop') %>% arrange(desc(NodeScore))
    fpr <- c(0,cumsum(grouped$fp)/n_neg)
    tpr <- c(0,cumsum(grouped$tp)/n_pos)
    auc <- sum(diff(fpr)*(head(tpr,-1)+tail(tpr,-1))/2)
    roc_rows[[length(roc_rows)+1]] <- tibble(seed=seed_id,comparison=comparison,fpr=fpr,tpr=tpr)
    auc_rows[[length(auc_rows)+1]] <- tibble(seed=seed_id,comparison=comparison,
                                             AUROC=auc,n_positive=n_pos,n_negative=n_neg)
    # Interpolate nonvertical ROC segments; max TPR at duplicate FPR knots.
    knots <- tibble(fpr=fpr,tpr=tpr) %>% group_by(fpr) %>% summarise(tpr=max(tpr),.groups='drop')
    roc_grid_rows[[length(roc_grid_rows)+1]] <- tibble(seed=seed_id,comparison=comparison,
                                                       fpr=fpr_grid,tpr=approx(knots$fpr,knots$tpr,xout=fpr_grid,rule=2)$y)
  }
}
auc_per_seed <- bind_rows(auc_rows)
auc_summary <- auc_per_seed %>% group_by(comparison) %>%
  summarise(n_datasets=n(),mean_AUROC=mean(AUROC),sd_AUROC=sd(AUROC),.groups='drop')
roc_mean <- bind_rows(roc_grid_rows) %>% group_by(comparison,fpr) %>%
  summarise(mean_tpr=mean(tpr),.groups='drop')
roc_mean <- bind_rows(tibble(comparison=names(roc_comparisons),fpr=0,mean_tpr=0),roc_mean) %>% arrange(comparison,fpr,mean_tpr)
roc_colors <- c('A vs B'='#E8734E','C vs B'='#704B95','A vs C'='#227C9D')
roc_labels <- setNames(paste0(auc_summary$comparison,'  AUROC: ',sprintf('%.2f',auc_summary$mean_AUROC)),auc_summary$comparison)
for(which_plot in c('All_Comparisons')) {
  comparisons <- if(which_plot=='A_vs_C') 'A vs C' else if(which_plot=='All_Comparisons') c('A vs B','C vs B','A vs C') else c('A vs B','C vs B')
  p_roc <- ggplot(roc_mean %>% filter(comparison %in% comparisons),aes(fpr,mean_tpr,color=comparison))+
    geom_abline(intercept=0,slope=1,linetype='dashed',color='black',linewidth=.7)+
    geom_line(linewidth=1.25)+
    scale_color_manual(values=roc_colors,breaks=comparisons,labels=roc_labels[comparisons])+
    scale_x_continuous(breaks=seq(0,1,.25),limits=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,.25),limits=c(0,1))+
    coord_equal()+
    labs(title=if(which_plot=='A_vs_C') 'Node-level AUROC: Set A versus Set C' else 'Node-level AUROC comparisons',
         subtitle='Cohort node scores; curves averaged across datasets',
         x='False positive rate',y='True positive rate',color=NULL)+
    figure_theme+theme(legend.position=c(.64,.15),legend.background=element_rect(fill='white',color=NA),
                       legend.text=element_text(size=11))
  save_plot(p_roc,paste0('Fig2_Node_ROC_',which_plot),7,6.2)
}
write_csv(bind_rows(roc_rows),file.path(plot_dir,'Table_Node_ROC_Per_Seed.csv'))
write_csv(roc_mean,file.path(plot_dir,'Table_Node_ROC_Mean.csv'))
write_csv(auc_per_seed,file.path(plot_dir,'Table_Node_AUROC_Per_Seed.csv'))
write_csv(auc_summary,file.path(plot_dir,'Table_Node_AUROC_Summary.csv'))

# =============================================================================
# Figure-2: top-node recovery across independent simulated datasets
# =============================================================================

  # Top-node counts, expected counts and class-size-adjusted recall.
  top_nodes <- nodes %>%
    group_by(seed) %>%
    arrange(median_patient_rank, node_idx, .by_group = TRUE) %>%
    slice_head(n = cohort_top_k) %>%
    ungroup()
  
  totals <- nodes %>%
    count(seed, set_membership, name = 'set_size')
  
  top_counts <- top_nodes %>%
    count(seed, set_membership, name = 'selected')
  
  top_counts <- totals %>%
    left_join(top_counts, by = c('seed', 'set_membership')) %>%
    mutate(selected = replace_na(selected, 0L)) %>%
    group_by(seed) %>%
    mutate(
      k = sum(selected),
      expected = k * set_size / sum(set_size),
      recall = selected / set_size,
      fold_enrichment = selected / expected
    ) %>%
    ungroup()
  
  write_csv(top_counts, file.path(output_dir, 'top_node_recovery.csv'))
  
  count_summary <- top_counts %>%
    group_by(set_membership) %>%
    summarise(
      mean_selected = mean(selected),
      expected = mean(expected),
      .groups = 'drop'
    )
  
  p_top <- ggplot(count_summary,
    aes(set_membership, mean_selected, fill = set_membership) ) +
    geom_col(alpha = 0.65) +
    geom_point(aes(y = expected), shape = 95, size = 14) +
    geom_point(
      data = top_counts,
      aes(y = selected), size = 2.5,
      position = position_jitter(width = 0.28, height = 0)
    ) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = paste0(
        'Top-', min(cohort_top_k, nrow(example_nodes)), ' cohort nodes'
      ),
      x = 'Node set',
      y = 'Number selected',
      caption = paste(
        'Dots = independent datasets;',
        'black marks = random-selection expectation.'
      )
    ) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none') 
    
  
  save_plot(p_top, 'Top_nodes_vs_expectation', width = 3.6, height = 4)
  


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
  # 1. Update factor levels and labels to match the x-axis text in the image
  example_nodes <- example_nodes %>%
    mutate(node_class = factor(set_membership, 
                               levels = c('A', 'C', 'B'),
                               labels = c('Exposure-scaled\n(Set-A)', 
                                          'Cancer-only\n(Set-C)', 
                                          'Background\n(Set-B)')))
  
  # 2. Define pairwise comparisons for the p-value brackets
  my_comparisons <- list(
    c('Exposure-scaled\n(Set-A)', 'Cancer-only\n(Set-C)'),
    c('Exposure-scaled\n(Set-A)', 'Background\n(Set-B)'),
    c('Cancer-only\n(Set-C)', 'Background\n(Set-B)')
  )
  
  # 3. Define the specific colors matching the image
  class_colors <- c('Exposure-scaled\n(Set-A)' = '#E79A84', # Peach/Coral
                    'Cancer-only\n(Set-C)' = '#9B82B5',     # Purple
                    'Background\n(Set-B)' = '#7A7A7A')      # Grey
  
  # 4. Build the plot using `normalized_priority` (the actual column in your data)
  p_node <- ggplot(example_nodes, aes(x = node_class, y = normalized_priority, fill = node_class)) +
    geom_violin(alpha = 0.1, trim = FALSE, color = "black") + 
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.05, size = 1.2, alpha = 0.6, color = "black") + 
    scale_fill_manual(values = class_colors) +
    labs(title = 'Cohort node-priority score by node class', 
         x = NULL, 
         y = 'Cohort node-priority score') + 
    theme_bw(base_size = 14) + # Light grid background
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 18),
          axis.text.x = element_text(size = 16, color = "grey30"),
          axis.title.y = element_text(size = 16)) +
    # Add Pairwise p-values
    stat_compare_means(comparisons = my_comparisons, 
                       method = "wilcox.test", 
                       label = "p.format",
                       tip.length = 0.01) +
    # Add Global Kruskal-Wallis p-value (dynamically adjusting y position)
    stat_compare_means(method = "kruskal.test", 
                       label.x = 0.6, 
                       label.y = max(example_nodes$normalized_priority, na.rm = TRUE) * 1.3)
  
  # Save plot 
  save_plot(p_node, 'Cohort_node_priority_by_class', width = 6.5, height = 7)
  
} else {
  p_node <- example_nodes %>% arrange(median_patient_rank, node_idx) %>% slice_head(n = 30) %>%
    ggplot(aes(reorder(node_label, normalized_priority), normalized_priority)) + geom_col(fill = '#337A96') +
    coord_flip() + labs(title = 'Top cohort nodes', x = NULL, y = 'Normalized rank priority')
  save_plot(p_node, 'Top_cohort_nodes', width = 8, height = 5)
}


# =============================================================================
# Supplement 3: cancer-transition node score (paired embedding change)
# =============================================================================
if (has_truth) {
  example_nodes <- example_nodes %>%
    mutate(node_class = factor(set_membership, 
                               levels = c('A', 'C', 'B'),
                               labels = c('Exposure-scaled\n(Set-A)', 
                                          'Cancer-only\n(Set-C)', 
                                          'Background\n(Set-B)')))
  
  my_comparisons <- list(
    c('Exposure-scaled\n(Set-A)', 'Cancer-only\n(Set-C)'),
    c('Exposure-scaled\n(Set-A)', 'Background\n(Set-B)'),
    c('Cancer-only\n(Set-C)', 'Background\n(Set-B)')
  )
  
  class_colors <- c('Exposure-scaled\n(Set-A)' = '#E79A84', 
                    'Cancer-only\n(Set-C)' = '#9B82B5',     
                    'Background\n(Set-B)' = '#7A7A7A')      
  
  p_cancer <- ggplot(example_nodes, aes(x = node_class, y = median_patient_transition_score, fill = node_class)) +
    geom_violin(alpha = 0.1, trim = FALSE, color = "black") + 
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.05, size = 1.2, alpha = 0.6, color = "black") + 
    scale_fill_manual(values = class_colors) +
    labs(title = 'Node Delta by simulated node class', 
         x = NULL, 
         y = 'Paired-transition magnitude') + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 18),
          axis.text.x = element_text(size = 16, color = "grey30"),
          axis.title.y = element_text(size = 16)) +
    stat_compare_means(comparisons = my_comparisons, 
                       method = "wilcox.test", 
                       label = "p.format",
                       tip.length = 0.01) +
    stat_compare_means(method = "kruskal.test", 
                       label.x = 0.6, 
                       label.y = max(example_nodes$median_patient_transition_score, na.rm = TRUE) * 1.2)
  
  save_plot(p_cancer, 'Supp3_node_transition_score', width = 8.5, height = 7)
  
} else {
  p_cancer <- ggplot(example_nodes, aes(median_patient_transition_score)) + 
    geom_histogram(bins = 40, fill = '#337A96', color = 'white') +
    labs(title = 'Paired-state node perturbation', 
         x = 'Median patient node-transition score', 
         y = 'Number of nodes')
  save_plot(p_cancer, 'Supp3_node_transition_score', width = 6.5, height = 7)
}


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
# Directed incoming/outgoing counts stay separate. These are model message channels, not measured interactions.
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
# Supplementary Figure: repetition controls and within-set recall
# =============================================================================

# Folder containing repetition_AUROC.csv and topk_recovery.csv
figure_dir <- file.path(dirname(manifest_file), 'figures')  # written by plot_simulation.py

repetitions <- read_csv(
  file.path(figure_dir, "repetition_AUROC.csv"),
  show_col_types = FALSE
)

recovery <- read_csv(
  file.path(figure_dir, "topk_recovery.csv"),
  show_col_types = FALSE
)

set.seed(42)  # Reproducible horizontal jitter only

# =============================================================================
# Panel A: A-versus-C recovery under three repetition controls
# =============================================================================

# First average repetitions within each independent simulated dataset.
# Each plotted dot therefore represents one dataset, not one model fit.
dataset_auc <- repetitions %>%
  group_by(seed, kind) %>%
  summarise(
    AUROC = mean(A_vs_C_AUROC),
    .groups = "drop"
  ) %>%
  mutate(
    kind = factor(
      kind,
      levels = c("initialization", "subset", "permutation")
    )
  )

auc_summary <- dataset_auc %>%
  group_by(kind) %>%
  summarise(
    datasets = n(),
    mean_auc = mean(AUROC),
    sd_auc = sd(AUROC),
    .groups = "drop"
  )

p_controls <- ggplot(dataset_auc, aes(kind, AUROC)) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.5,
    color = "#337A96"
  ) +
  geom_errorbar(
    data = auc_summary,
    aes(x = kind, ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_point(
    data = auc_summary,
    aes(x = kind, y = mean_auc),
    inherit.aes = FALSE,
    shape = 95,
    size = 9
  ) +
  scale_x_discrete(
    labels = c(
      initialization = "Different\ninitializations",
      subset = "Patient\nsubsets",
      permutation = "Permuted\nexposure"
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Three repetition controls",
    x = NULL,
    y = "Mean single-run A-versus-C AUROC"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# =============================================================================
# Panel B: fraction of each node set recovered in the top-50 list
# =============================================================================

# recall = selected nodes from a set / total nodes in that set
recall_data <- recovery %>%
  mutate(set = factor(set, levels = c("A", "C", "B")))

recall_summary <- recall_data %>%
  group_by(set) %>%
  summarise(
    datasets = n(),
    mean_recall = mean(recall),
    sd_recall = sd(recall),
    .groups = "drop"
  )

p_recall <- ggplot(recall_data, aes(set, recall, color = set)) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.5
  ) +
  geom_errorbar(
    data = recall_summary,
    aes(
      x = set,
      ymin = mean_recall - sd_recall,
      ymax = mean_recall + sd_recall
    ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_point(
    data = recall_summary,
    aes(x = set, y = mean_recall),
    inherit.aes = FALSE,
    shape = 95,
    size = 9
  ) +
  scale_color_manual(
    values = c(A = "#E8734E", C = "#704B95", B = "#8A8A8A")
  ) +
  scale_x_discrete(
    labels = c(A = "Set A", C = "Set C", B = "Set B")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Recovery adjusted for set size",
    x = "Node set",
    y = "Top-50 recall within each set"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# =============================================================================
# Combine panels and export
# =============================================================================

supp5 <- (p_controls | p_recall) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Dots represent independent simulated datasets;",
      "horizontal marks and error bars show mean and SD across datasets.",
      "\nPanel A averages model repetitions within each dataset first.",
      "Permutation results are descriptive."
    )
  )

ggsave(file.path(figure_dir, "Supplement_5_robustness.pdf"),supp5, width = 11, height = 5.5)




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

# ============================================================
# Model-defined edge attention changes: new exports, original plot theme
# ============================================================
# manifest_file is set from the command line (argument 4)
plot_dir <- output_dir
association_mode <- 'signed_change' # Primary: signed edge-weight change vs exposure.
# 'magnitude' reproduces the older question: absolute edge change vs exposure.
# Select the definition BEFORE inspecting which gives the preferred A/C separation.
top_fraction <- 0.10
histogram_bins <- 40
minimum_change_sd <- 1e-12 # Numerical nonvariation filter, not significance.
# Average reciprocal directions into one connection BEFORE association calculation.
# Set FALSE only if directed channels are the intended edge unit for every plot.
undirected <- TRUE
show_seed_dots <- TRUE

args <- commandArgs(trailingOnly=TRUE)
# arguments are parsed at the top of the script
# association_mode is fixed above (signed_change)
stopifnot(association_mode %in% c('signed_change','magnitude'))
dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
manifest <- read_csv(manifest_file,show_col_types=FALSE)
stopifnot(all(c('seed','data_dir','analysis_dir') %in% names(manifest)))
if(anyDuplicated(manifest$seed)) stop('Use one row per independent dataset.')
class_order <- c('A','C','B')
class_colors <- c(A='#E8734E',C='#704B95',B='#8A8A8A')
class_labels <- c(A='Exposure-associated\n(Set-A)',C='Cancer only\n(Set-C)',B='Background\n(Set-B)')
set.seed(42) # Plot jitter only.

# ============================================================
# Small plotting and correlation helpers
# ============================================================
safe_spearman <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<4 || sd(x[ok])<=minimum_change_sd || sd(y[ok])<=minimum_change_sd) return(NA_real_)
  suppressWarnings(cor(x[ok],y[ok],method='spearman'))
}
save_plot <- function(p,name,width=7,height=5) {
  ggsave(file.path(plot_dir,paste0(name,'.pdf')),p,width=width,height=height)
  ggsave(file.path(plot_dir,paste0(name,'.png')),p,width=width,height=height,dpi=300)
}
figure_theme <- theme_bw(base_size=13)+theme(plot.title=element_text(face='bold'),
                                             axis.text.x=element_text(size=11),legend.position='none')
edge_list <- list();count_list <- list();degree_list <- list();hist_list <- list();selection_list <- list()

# ============================================================
# Compute patient-specific edge change and cohort exposure association
# ============================================================
for(i in seq_len(nrow(manifest))) {
  id <- manifest$seed[i]; data_dir <- manifest$data_dir[i]; analysis_dir <- manifest$analysis_dir[i]
  config <- fromJSON(file.path(analysis_dir,'config.json'))
  truth <- read_csv(file.path(data_dir,'sim_ground_truth_nodes.csv'),show_col_types=FALSE) %>%
    select(node_idx,node_label,set_membership)
  if(anyDuplicated(truth$node_idx) || any(!truth$set_membership %in% class_order)) stop('Expected current node labels A/C/B and unique node indices.')
  exposure <- read_csv(file.path(data_dir,'metadata_samples.csv'),show_col_types=FALSE) %>%
    transmute(patient_id=as.character(patient_id),exposure=.data[[config$exposure_col]]) %>% distinct()
  if(anyDuplicated(exposure$patient_id)) stop('Inconsistent paired exposure values.')
  
  # This is the median across initializations of signed gated-attention changes.
  # Do not substitute median(tumor attention)-median(normal attention): not equivalent.
  ep <- read_csv(file.path(analysis_dir,'patient_edge_scores.csv.gz'),show_col_types=FALSE,
                 progress=FALSE,col_types=cols_only(patient_id=col_character(),edge_id=col_integer(),
                                                    source_idx=col_integer(),target_idx=col_integer(),median_gated_attention_change=col_double()))
  if(anyDuplicated(ep[c('patient_id','edge_id')])) stop('Duplicate patient/edge rows.')
  if(!all(ep$patient_id %in% exposure$patient_id)) stop('Some scored patients lack exposure metadata.')
  ep <- ep %>% filter(source_idx!=target_idx) %>%
    mutate(u=if(undirected) pmin(source_idx,target_idx) else source_idx,
           v=if(undirected) pmax(source_idx,target_idx) else target_idx) %>%
    group_by(patient_id,u,v) %>% summarise(delta_weight=mean(median_gated_attention_change),.groups='drop') %>%
    left_join(exposure,by='patient_id') %>% mutate(attention_change_magnitude=abs(delta_weight))
  if(any(!is.finite(ep$delta_weight)) || any(!is.finite(ep$exposure))) stop('Nonfinite edge changes or exposures; fix input provenance instead of silently dropping rows.')
  expected_patients <- n_distinct(ep$patient_id)
  if(any((ep %>% count(u,v))$n!=expected_patients)) stop('Edges do not have the same scored patient coverage.')
  write_csv(ep,file.path(plot_dir,paste0('Patient_Edge_AttentionChange_seed_',id,'.csv.gz')))
  edge <- ep %>% group_by(u,v) %>% summarise(
    n_patients=n(),mean_attention_change_magnitude=mean(attention_change_magnitude),
    median_attention_change_magnitude=median(attention_change_magnitude),
    rho_signed_change=safe_spearman(delta_weight,exposure),
    rho_magnitude=safe_spearman(attention_change_magnitude,exposure),.groups='drop') %>%
    mutate(ExposureEdgeScore=if(association_mode=='signed_change') rho_signed_change else rho_magnitude,
           EdgePriorityScore=abs(ExposureEdgeScore),seed=id) %>%
    left_join(truth %>% transmute(u=node_idx,src_set=set_membership),by='u') %>%
    left_join(truth %>% transmute(v=node_idx,dst_set=set_membership),by='v') %>%
    mutate(EdgeClass=if_else(src_set==dst_set,src_set,'Mixed'))
  if(anyNA(edge$src_set) || anyNA(edge$dst_set)) stop('Unmapped edge endpoint; never assume missing truth is background.')
  
  # ============================================================
  # Select top edges globally; then report within-set and mixed edges
  # ============================================================
  requested_k <- max(1,ceiling(top_fraction*nrow(edge)))
  selected <- edge %>% filter(is.finite(EdgePriorityScore),EdgePriorityScore>0) %>%
    arrange(desc(EdgePriorityScore),u,v) %>% slice_head(n=requested_k) %>% select(u,v) %>% mutate(selected=TRUE)
  edge <- edge %>% left_join(selected,by=c('u','v')) %>% mutate(selected=replace_na(selected,FALSE))
  edge_list[[i]] <- edge
  selection_list[[i]] <- tibble(seed=id,total_edges=nrow(edge),requested_k=requested_k,
                                selected_edges=sum(edge$selected),selected_within_set=sum(edge$selected & edge$EdgeClass!='Mixed'),
                                selected_mixed=sum(edge$selected & edge$EdgeClass=='Mixed'),
                                undefined_correlations=sum(!is.finite(edge$ExposureEdgeScore)))
  eligible_total <- sum(is.finite(edge$EdgePriorityScore)&edge$EdgePriorityScore>0)
  counts <- edge %>% group_by(EdgeClass) %>% summarise(available_edges=n(),n_top_edges=sum(selected),
                                                       eligible_edges=sum(is.finite(EdgePriorityScore)&EdgePriorityScore>0),.groups='drop') %>%
    complete(EdgeClass=c(class_order,'Mixed'),fill=list(available_edges=0,n_top_edges=0,eligible_edges=0)) %>%
    mutate(seed=id,
           expected_all_graph=sum(n_top_edges)*available_edges/sum(available_edges),
           expected_eligible=if(eligible_total>0) sum(n_top_edges)*eligible_edges/eligible_total else 0)
  within_selected <- sum(counts$n_top_edges[counts$EdgeClass!='Mixed'])
  counts <- counts %>% mutate(percent_of_all_selected=if(sum(n_top_edges)>0) 100*n_top_edges/sum(n_top_edges) else NA_real_,
                              percent_of_within_selected=if_else(EdgeClass!='Mixed' & within_selected>0,100*n_top_edges/within_selected,NA_real_))
  count_list[[i]] <- counts
  
  # ============================================================
  # Node selected degree and weighted degree
  # ============================================================
  # Original-style panel uses within-class edges. Also save all-edge sensitivity.
  for(scope_name in c('within_class','all_edges')) {
    chosen <- edge %>% filter(selected)
    if(scope_name=='within_class') chosen <- chosen %>% filter(EdgeClass!='Mixed')
    incident <- bind_rows(chosen %>% transmute(node_idx=u,edge_weight=EdgePriorityScore),
                          chosen %>% transmute(node_idx=v,edge_weight=EdgePriorityScore)) %>%
      group_by(node_idx) %>% summarise(SelectedDegree=n(),WeightedDegree=sum(edge_weight),.groups='drop')
    nd <- truth %>% left_join(incident,by='node_idx') %>%
      mutate(SelectedDegree=replace_na(SelectedDegree,0L),WeightedDegree=replace_na(WeightedDegree,0),
             seed=id,scope=scope_name)
    degree_list[[length(degree_list)+1]] <- nd
  }
  
  # ============================================================
  # Histogram bins: actual signed correlations, not rank percentiles
  # ============================================================
  hist <- edge %>% filter(EdgeClass!='Mixed',is.finite(ExposureEdgeScore)) %>%
    mutate(bin=pmin(histogram_bins-1L,pmax(0L,floor((ExposureEdgeScore+1)/2*histogram_bins)))) %>%
    count(EdgeClass,bin,name='n_edges') %>%
    complete(EdgeClass=class_order,bin=0:(histogram_bins-1),fill=list(n_edges=0)) %>% mutate(seed=id)
  hist_list[[i]] <- hist
}
all_edges <- bind_rows(edge_list);counts <- bind_rows(count_list)
node_degrees <- bind_rows(degree_list);histograms <- bind_rows(hist_list)
write_csv(all_edges,file.path(plot_dir,'Table_Edges_With_AttentionChange_And_EdgeScore.csv'))
write_csv(counts,file.path(plot_dir,'Table_Top10_Edge_Composition.csv'))
write_csv(node_degrees,file.path(plot_dir,'Table_Node_WeightedDegree_By_TruthClass.csv'))
write_csv(bind_rows(selection_list),file.path(plot_dir,'Table_Selected_And_Mixed_Edge_Counts.csv'))

# ============================================================
# A. Composition of top top-scored edges: original three bars
# ============================================================
composition <- counts %>% filter(EdgeClass!='Mixed') %>% group_by(EdgeClass) %>% summarise(
  sd_count=sd(n_top_edges),mean_count=mean(n_top_edges),
  mean_percent=if(all(is.na(percent_of_within_selected))) NA_real_ else mean(percent_of_within_selected,na.rm=TRUE),
  mean_expected=mean(expected_all_graph),.groups='drop')
p_top <- ggplot(composition,aes(factor(EdgeClass,levels=class_order),mean_count,fill=EdgeClass)) +
  geom_col(width=.68,alpha=.85)+
  geom_errorbar(aes(ymin=mean_count-sd_count,ymax=mean_count+sd_count),width=.18,linewidth=.7,na.rm=TRUE)+
  geom_point(aes(y=mean_expected),shape=95,size=9)+
  geom_text(aes(y=mean_count+coalesce(sd_count,0),label=paste0('mean n=',round(mean_count,1),'\n(',round(mean_percent,1),'%)')),
            vjust=-.35,size=3.8)+
  scale_fill_manual(values=class_colors)+scale_x_discrete(labels=class_labels)+
  scale_y_continuous(expand=expansion(mult=c(.02,.25)))+
  labs(title=paste0('Top ',round(100*top_fraction),'% top-scored edges'),
       subtitle='Within-set edges; ranked by absolute exposure association',x=NULL,y='Mean number of selected edges',
       caption='Percentages: within-set selected edges only. Mixed edges reported separately.\nBlack marks: expected counts under uniform edge selection.')+figure_theme
if(show_seed_dots) p_top <- p_top+geom_point(data=counts %>% filter(EdgeClass!='Mixed'),aes(y=n_top_edges),
                                             position=position_jitter(width=.055,height=0),size=1.3,alpha=.65)
save_plot(p_top,'Fig_A_Top10_Edge_Composition',8,5.6)

# ============================================================
# Stacked counts: class composition of the top 10% selected edges
# ============================================================
# Only selected within-class edges: A-A, C-C, B-B. Mixed edges excluded.
# Use the EXISTING selection; do not rerank or select within each class.
stack_order <- c('A','C','B')
stack_colors <- class_colors
stack_labels <- c(A='Exposure scaled (A-A)',C='Cancer only (C-C)',B='Background (B-B)')
stack_counts <- counts %>% filter(EdgeClass!='Mixed') %>% group_by(EdgeClass) %>%
  summarise(mean_count=mean(n_top_edges),.groups='drop') %>%
  mutate(percent=if(sum(mean_count)>0) 100*mean_count/sum(mean_count) else NA_real_,
         label=if_else(mean_count>0,paste0(round(mean_count,1),'\n(',round(percent,1),'%)'),''))
p_stack <- ggplot(stack_counts,aes(x='Top 10% edges',y=mean_count,
                                   fill=factor(EdgeClass,levels=stack_order)))+
  geom_col(width=.5,color='white',linewidth=.4)+
  geom_text(aes(label=label),position=position_stack(vjust=.5),size=3.8)+
  scale_fill_manual(values=stack_colors,breaks=stack_order,labels=stack_labels)+
  scale_y_continuous(expand=expansion(mult=c(0,.06)))+
  labs(title='Within-class composition of top 10% ranked edges',x=NULL,
       y='Mean number of selected edges per dataset',fill='Edge class',
       caption='Selected globally using the existing edge priority score.\nOnly A-A, C-C, B-B edges; mixed/boundary edges excluded.')+
  figure_theme+theme(legend.position='right')
save_plot(p_stack,'Fig_Top10_Edge_Counts_Stacked',7,5.5)
write_csv(stack_counts,file.path(plot_dir,'Table_Top10_Edge_Counts_Stacked.csv'))
# Individual datasets: counts rather than pooled edges or just one seed.
p_stack_seed <- ggplot(counts %>% filter(EdgeClass!='Mixed'),aes(x=factor(seed),y=n_top_edges,
                                                                 fill=factor(EdgeClass,levels=stack_order)))+
  geom_col(width=.7,color='white',linewidth=.25)+
  scale_fill_manual(values=stack_colors,breaks=stack_order,labels=stack_labels)+
  scale_y_continuous(expand=expansion(mult=c(0,.06)))+
  labs(title='Within-class composition of top 10% ranked edges',x='Simulation seed',
       y='Number of selected edges',fill='Edge class')+
  figure_theme+theme(legend.position='right')
save_plot(p_stack_seed,'Fig_Top10_Edge_Counts_Stacked_Per_Seed',8,5.5)

# ============================================================
# B. Distribution of exposure association of edge attention changes
# ============================================================
hist_mean <- histograms %>% group_by(EdgeClass,bin) %>% summarise(mean_edges=mean(n_edges),.groups='drop') %>%
  mutate(score_midpoint=-1+(bin+.5)*2/histogram_bins)
p_hist <- ggplot(hist_mean,aes(score_midpoint,mean_edges,fill=factor(EdgeClass,levels=class_order)))+
  geom_col(width=2/histogram_bins,alpha=.85,color='black',linewidth=.3,position='identity')+
  scale_fill_manual(values=class_colors,labels=c('Exposure-associated','Cancer only','Background'))+
  scale_x_continuous(limits=c(-1,1),breaks=seq(-1,1,.5))+
  labs(title='Exposure association of edge attention changes',
       subtitle=if(association_mode=='signed_change') 'Signed edge-weight change vs exposure' else 'Absolute edge-weight change vs exposure',
       x='Exposure Edge Score (Spearman correlation)',y='Mean number of edges per dataset',fill='Edge class')+
  figure_theme+theme(legend.position='right')
save_plot(p_hist,'Fig_B_ExposureEdgeScore_Distribution',8,5.2)
save_plot((p_top/p_hist)+plot_annotation(tag_levels='A'),'Supp_Figure_4_Edge_Attention_Change',9,11)

# ============================================================
# B2. Older magnitude-based exposure association (additional plot)
# ============================================================
# Correlate ABSOLUTE paired change with exposure; do NOT take abs(correlation).
# Keep the primary selection and degree definitions unchanged.
magnitude_bins <- all_edges %>%
  filter(EdgeClass!='Mixed',is.finite(rho_magnitude)) %>%
  mutate(bin=pmin(histogram_bins-1L,pmax(0L,floor((rho_magnitude+1)/2*histogram_bins)))) %>%
  count(seed,EdgeClass,bin,name='n_edges') %>%
  complete(seed=manifest$seed,EdgeClass=class_order,bin=0:(histogram_bins-1),fill=list(n_edges=0))
magnitude_hist <- magnitude_bins %>% group_by(EdgeClass,bin) %>%
  summarise(mean_edges=mean(n_edges),.groups='drop') %>%
  mutate(score_midpoint=-1+(bin+.5)*2/histogram_bins)
p_magnitude <- ggplot(magnitude_hist,
                      aes(score_midpoint,mean_edges,fill=factor(EdgeClass,levels=class_order)))+
  geom_col(width=2/histogram_bins,alpha=.85,color='black',linewidth=.3,position='identity')+
  scale_fill_manual(values=class_colors,labels=c('Exposure-associated','Cancer only','Background'))+
  scale_x_continuous(limits=c(-1,1),breaks=seq(-1,1,.5))+
  labs(title='Exposure association of edge attention changes',
       subtitle='Magnitude of paired edge-weight change vs exposure',
       x='Exposure Edge Score (magnitude-based Spearman correlation)',
       y='Mean number of edges per dataset',fill='Edge class')+
  figure_theme+theme(legend.position='right')
save_plot(p_magnitude,'Fig_B2_Magnitude_ExposureEdgeScore_Distribution',8,5.2)
write_csv(magnitude_bins,file.path(plot_dir,'Table_Magnitude_Histogram_Per_Seed.csv'))
write_csv(magnitude_hist,file.path(plot_dir,'Table_Magnitude_Histogram_Mean_Counts.csv'))

# ============================================================
# All within-class edges: overlapping edge-score densities
# ============================================================
# Edge score = absolute exposure correlation (existing EdgePriorityScore).
# No top-10% filtering. Mixed/boundary edges excluded.
# Estimate each seed/class separately, then average curves across seeds.
# Reflection at 0 and 1 corrects boundary leakage for this bounded score.
density_bandwidth <- 0.05 # Same bandwidth for all classes and datasets.
density_grid <- seq(0,1,length.out=401)
density_input <- all_edges %>% filter(EdgeClass!='Mixed',is.finite(EdgePriorityScore))
density_curves <- list()
for(seed_id in manifest$seed) {
  for(class_id in class_order) {
    values <- density_input %>% filter(seed==seed_id,EdgeClass==class_id) %>% pull(EdgePriorityScore)
    if(length(values)==0) stop('Cannot estimate density: no finite edges for seed ',seed_id,' class ',class_id)
    curve <- vapply(density_grid,function(x) mean(
      dnorm((x-values)/density_bandwidth)+dnorm((x+values)/density_bandwidth)+
        dnorm((x-(2-values))/density_bandwidth))/density_bandwidth,numeric(1))
    # Normalize numerical area to one for every seed/class.
    area <- sum((head(curve,-1)+tail(curve,-1))/2*diff(density_grid))
    density_curves[[length(density_curves)+1]] <- tibble(seed=seed_id,EdgeClass=class_id,
                                                         edge_score=density_grid,density=curve/area,n_edges=length(values))
  }
}
density_per_seed <- bind_rows(density_curves)
density_mean <- density_per_seed %>% group_by(EdgeClass,edge_score) %>%
  summarise(density=mean(density),.groups='drop')
p_density <- ggplot(density_mean,aes(edge_score,density,
                                     fill=factor(EdgeClass,levels=class_order),color=factor(EdgeClass,levels=class_order)))+
  geom_area(position='identity',alpha=.28,linewidth=.65)+
  scale_fill_manual(values=class_colors,labels=c('Set A (A-A)','Set C (C-C)','Set B (B-B)'))+
  scale_color_manual(values=class_colors,labels=c('Set A (A-A)','Set C (C-C)','Set B (B-B)'))+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=expansion(mult=c(0,0)))+
  scale_y_continuous(expand=expansion(mult=c(0,.05)))+
  labs(title='Exposure association of edge attention changes',
       subtitle='All within-class edges; no top-10% selection',
       x='Edge score (absolute exposure correlation)',y='Density',fill='Edge class',color='Edge class',
       caption='Each class has unit area; density curves averaged equally across simulation datasets.\nMixed/boundary edges excluded. Common bandwidth = 0.05; boundary correction at 0 and 1.')+
  figure_theme+theme(legend.position='right')
save_plot(p_density,'Fig_All_WithinClass_EdgeScore_Density',8,5.5)
write_csv(density_input %>% select(seed,u,v,EdgeClass,EdgePriorityScore,ExposureEdgeScore),
          file.path(plot_dir,'Table_All_WithinClass_EdgeScores_For_Density.csv'))
write_csv(density_mean,file.path(plot_dir,'Table_EdgeScore_Density_Mean.csv'))
write_csv(density_per_seed,file.path(plot_dir,'Table_EdgeScore_Density_Per_Seed.csv'))

# ============================================================
# C. Average weighted degree: original style
# ============================================================
# Average within each dataset first; use between-dataset SD, NOT node-wise SEM.
seed_degree <- node_degrees %>% group_by(seed,scope,set_membership) %>% summarise(
  mean_weighted_degree=mean(WeightedDegree),mean_selected_degree=mean(SelectedDegree),.groups='drop')
degree_summary <- seed_degree %>% group_by(scope,set_membership) %>% summarise(
  sd_weighted=sd(mean_weighted_degree),sd_degree=sd(mean_selected_degree),
  mean_weighted=mean(mean_weighted_degree),mean_degree=mean(mean_selected_degree),.groups='drop')
write_csv(seed_degree,file.path(plot_dir,'Table_Node_Degree_Per_Seed.csv'))
write_csv(degree_summary,file.path(plot_dir,'Table_Node_WeightedDegree_Summary.csv'))
write_csv(composition,file.path(plot_dir,'Table_Edge_Composition_Summary.csv'))
write_csv(hist_mean,file.path(plot_dir,'Table_Histogram_Mean_Counts.csv'))
for(scope_name in c('within_class','all_edges')) {
  d <- degree_summary %>% filter(scope==scope_name)
  p_degree <- ggplot(d,aes(factor(set_membership,levels=class_order),mean_weighted,fill=set_membership))+
    geom_col(width=.68,alpha=.85)+
    geom_errorbar(aes(ymin=mean_weighted-sd_weighted,ymax=mean_weighted+sd_weighted),width=.18,linewidth=.8,na.rm=TRUE)+
    geom_text(aes(y=mean_weighted+coalesce(sd_weighted,0),label=round(mean_weighted,2)),vjust=-.45,size=4)+
    scale_fill_manual(values=class_colors)+scale_x_discrete(labels=class_labels)+
    scale_y_continuous(expand=expansion(mult=c(.02,.2)))+
    labs(title='Weighted degree among top 10% exposure-associated edges by simulated node class',
         subtitle=if(scope_name=='within_class') 'Selected within-set edges only' else 'All selected edges, including mixed edges',
         x=NULL,y='Average weighted degree',
         caption='Edge weight = absolute exposure correlation.\nBars: dataset means; error bars: between-dataset SD.')+figure_theme
  if(show_seed_dots) p_degree <- p_degree+geom_point(data=seed_degree %>% filter(scope==scope_name),
                                                     aes(y=mean_weighted_degree),position=position_jitter(width=.055,height=0),size=1.3,alpha=.65)
  save_plot(p_degree,paste0('Supp_Figure_5_Weighted_Degree_',scope_name),9,5.8)
}


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
  'Attention and message changes on a fixed topology are not measured changes in biological interactions.',
  'No CpG-to-gene relationship validation is performed.',
  'Per-epoch correlation and held-out validation are unavailable in these exports.',
  'Latent magnitude checks, when available, are read from simulator diagnostics.'
), file.path(output_dir, 'analysis_notes.txt'))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'R_session_info.txt'))
message('Saved plots and supporting tables to: ', normalizePath(output_dir))

