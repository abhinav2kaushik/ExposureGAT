# Rscript plot_label_controls.R LABEL_CONTROLS_DIR
# Figures for src/analysis/label_controls.py: ranking agreement (real vs shuffled labels),
# exposure-aligned gene counts and fixed top 5% agreement.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)})
a <- commandArgs(trailingOnly = TRUE); if (length(a) != 1) stop('Supply the label-controls folder')
out <- a[1]
stats <- read_csv(file.path(out, 'statistical_comparisons.csv'), show_col_types = FALSE)
p_label <- function(x) paste0('Holm-adjusted p = ', format(signif(x, 2), scientific = FALSE))
theme_set(theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold')))
save <- function(p, name, w, h) { ggsave(file.path(out, paste0(name, '.pdf')), p, width = w, height = h)
  ggsave(file.path(out, paste0(name, '.png')), p, width = w, height = h, dpi = 180) }
fmt <- function(x) x %>% pivot_longer(c(rank_spearman, top100_jaccard), names_to = 'metric', values_to = 'agreement') %>%
  mutate(comparison = factor(comparison, levels = c('Real vs real', 'Real vs shuffled')),
         metric = recode(metric, rank_spearman = 'Rank Spearman', top100_jaccard = 'Top-100 Jaccard'))
pat <- fmt(read_csv(file.path(out, 'patient_ranking_agreement.csv'), show_col_types = FALSE))
coh <- fmt(read_csv(file.path(out, 'cohort_ranking_agreement.csv'), show_col_types = FALSE))
col <- c('Real vs real' = '#388DAF', 'Real vs shuffled' = '#AF62A2')
ann <- stats %>% filter(level %in% c('patient', 'cohort')) %>%
  mutate(metric = recode(metric, rank_spearman = 'Rank Spearman', top100_jaccard = 'Top-100 Jaccard'), label = p_label(p_holm))
p1 <- ggplot(pat, aes(comparison, agreement)) + geom_line(aes(group = patient_id), color = 'grey75', alpha = .6) +
  geom_boxplot(aes(fill = comparison), alpha = .5, outlier.shape = NA, width = .5) +
  geom_point(aes(color = comparison), position = position_jitter(width = .07, height = 0, seed = 42), alpha = .7, size = 1.5) +
  facet_wrap(~metric, scales = 'free_y') + scale_fill_manual(values = col) + scale_color_manual(values = col) +
  geom_text(data = filter(ann, level == 'patient'), aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.3, size = 4) +
  labs(title = 'A. Patient-specific rankings (weighted node score x |signed contribution|)', x = NULL, y = 'Mean agreement over fit pairs') +
  theme(legend.position = 'none')
p2 <- ggplot(coh, aes(comparison, agreement)) + geom_boxplot(aes(fill = comparison), alpha = .5, outlier.shape = NA, width = .5) +
  geom_point(aes(color = comparison), position = position_jitter(width = .09, seed = 42), alpha = .5, size = 1.3) +
  facet_wrap(~metric, scales = 'free_y') + scale_fill_manual(values = col) + scale_color_manual(values = col) +
  geom_text(data = filter(ann, level == 'cohort'), aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.3, size = 4) +
  labs(title = 'B. Cohort weighted node-score rankings', x = NULL, y = 'Agreement per fit pair') + theme(legend.position = 'none')
save(p1 / p2, 'Ranking_agreement_real_vs_shuffled', 13, 11)

counts <- read_csv(file.path(out, 'alignment_counts.csv'), show_col_types = FALSE) %>%
  mutate(condition = factor(if_else(kind == 'initialization', 'Real labels', 'Shuffled labels'), levels = c('Real labels', 'Shuffled labels')))
pal <- c('Real labels' = '#388DAF', 'Shuffled labels' = '#AF62A2')
p3 <- ggplot(filter(counts, mode == 'Positive alignment', threshold == 0.3), aes(condition, count, fill = condition)) +
  geom_boxplot(width = .5, alpha = .5, outlier.shape = NA) + geom_point(position = position_jitter(width = .08, seed = 42), size = 2.5) +
  scale_fill_manual(values = pal) +
  annotate('text', x = 1.5, y = Inf, label = p_label(stats$p_holm[stats$metric == 'Positive alignment: 0.3']), vjust = 1.3, size = 4.5) +
  labs(title = "Genes aligned with each fit's own exposure labels (Spearman rho >= 0.30)", x = NULL, y = 'Number of genes') +
  theme(legend.position = 'none')
save(p3, 'Exposure_aligned_gene_counts', 10, 7)
counts %>% group_by(mode, threshold, condition) %>% summarise(mean = mean(count), median = median(count), sd = sd(count), .groups = 'drop') %>%
  write_csv(file.path(out, 'alignment_count_summary.csv'))

fx <- read_csv(file.path(out, 'fixed_top5_correlations.csv'), show_col_types = FALSE)
p4 <- ggplot(fx, aes(permutation, spearman)) + geom_hline(yintercept = 0, color = 'grey60', linetype = 'dashed') +
  geom_line(color = '#AF62A2', linewidth = .8) + geom_point(color = '#AF62A2', size = 2.8) +
  scale_x_continuous(breaks = fx$permutation) + coord_cartesian(ylim = c(-1, 1)) +
  labs(title = 'Top 5% node agreement', subtitle = paste0('Fixed top 5% genes of the real-label consensus (n = ', fx$genes[1],
       '); mean Spearman rho = ', sprintf('%.3f', mean(fx$spearman, na.rm = TRUE))), x = 'Exposure permutation', y = 'Spearman correlation coefficient')
save(p4, 'Fixed_top5_agreement_with_shuffled_fits', 12, 6.5)
