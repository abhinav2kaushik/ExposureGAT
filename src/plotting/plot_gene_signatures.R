# Rscript plot_gene_signatures.R GENE_SIGNATURE_DIR WEIGHTED_NODE_SCORE_DIR
# Bubble plots of signature representation (real vs shuffled labels) and weighted node-score
# densities of signature genes within the top 5/10/20% of real-label fits.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
a <- commandArgs(trailingOnly = TRUE); if (length(a) != 2) stop('Supply gene-signature and weighted-node-score folders')
out <- a[1]
s <- read_csv(file.path(out, 'signature_summary.csv'), show_col_types = FALSE)
background <- first(read_csv(file.path(out, 'signature_counts_by_fit.csv'), show_col_types = FALSE)$background)
colors <- c('Smoking/tobacco-related' = '#E87854', 'Smoking and lung cancer' = '#33869D', 'Lung-cancer-related' = '#76549A')
d <- s %>% filter(kind %in% c('initialization', 'permutation'), scored_genes > 0) %>%
  mutate(condition = factor(if_else(kind == 'initialization', 'Real labels', 'Shuffled labels'), levels = c('Real labels', 'Shuffled labels')),
         percent = 100 * mean_fraction, label = paste0(sprintf('%.1f', mean_selected), '/', scored_genes))
bubble <- function(z, name) {
  if (!nrow(z)) return(invisible(NULL))
  z$pathway <- factor(z$pathway, levels = rev(sort(unique(z$pathway))))
  p <- ggplot(z, aes(percent, pathway)) + geom_vline(xintercept = 100 * background, linetype = 'dashed', color = 'grey45') +
    geom_point(aes(size = mean_selected, color = category), alpha = .85) +
    geom_text(aes(label = label), nudge_x = 1.2, hjust = 0, size = 3.2) + facet_grid(. ~ condition) +
    scale_size_area(max_size = 8, name = 'Mean selected\ngenes per fit') + scale_color_manual(values = colors, name = 'Signature type') +
    scale_x_continuous(limits = c(0, max(100 * background, max(z$percent)) + 10)) +
    scale_y_discrete(labels = function(x) gsub('_', ' ', x)) + theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 9)) +
    labs(x = 'Mean selected genes (% of scored signature genes)', y = NULL)
  h <- max(5, .3 * length(unique(z$pathway)) + 2)
  ggsave(file.path(out, paste0(name, '.pdf')), p, width = 16, height = h, limitsize = FALSE)
  ggsave(file.path(out, paste0(name, '.png')), p, width = 16, height = h, dpi = 150, limitsize = FALSE)
}
bubble(filter(d, category != 'Lung-cancer-related'), 'Smoking_signatures')
bubble(filter(d, category == 'Lung-cancer-related'), 'Lung_cancer_signatures')

w <- read_csv(file.path(a[2], 'weighted_node_scores_by_fit.csv.gz'), show_col_types = FALSE, col_types = cols(entrez_id = 'c')) %>%
  filter(kind == 'initialization') %>% group_by(run) %>% arrange(desc(score), node_idx, .by_group = TRUE) %>%
  mutate(rank = row_number(), n_genes = n()) %>% ungroup()
sets <- read_csv(file.path(out, 'gene_set_membership.csv'), col_types = cols(.default = 'c')) %>% distinct(category, entrez_id)
sel <- bind_rows(lapply(c(.05, .10, .20), function(f) w %>% filter(rank <= ceiling(f * n_genes)) %>% mutate(fraction = f))) %>%
  inner_join(sets, by = 'entrez_id', relationship = 'many-to-many')
curves <- sel %>% group_by(category, run, fraction) %>% group_modify(~{
  v <- .x$score; dd <- density(c(v, -v, 2 - v), bw = .025, from = 0, to = 1, n = 512); tibble(score = dd$x, density = 3 * dd$y) }) %>% ungroup() %>%
  group_by(category, fraction, score) %>% summarise(density = mean(density), .groups = 'drop')
write_csv(curves, file.path(out, 'weighted_score_density_data.csv'))
pal <- c('Top 5%' = '#CF5D3A', 'Top 10%' = '#287D8E', 'Top 20%' = '#8263A7')
for (cat in unique(curves$category)) {
  x <- curves %>% filter(category == cat) %>% mutate(selection = factor(fraction, levels = c(.05, .10, .20), labels = names(pal)))
  p <- ggplot(x, aes(score, density, color = selection, fill = selection)) + geom_area(position = 'identity', alpha = .12) +
    geom_line(linewidth = 1.1) + scale_color_manual(values = pal) + scale_fill_manual(values = pal) + coord_cartesian(xlim = c(0, .6)) +
    labs(title = paste0('Distribution of weighted node score (', tolower(cat), ' genes)'), x = 'Weighted node score', y = 'Density', color = NULL, fill = NULL) +
    theme_bw(base_size = 15) + theme(panel.grid.minor = element_blank(), legend.position = 'right', plot.title = element_text(face = 'bold'))
  name <- paste0('Weighted_score_density_', gsub('[^A-Za-z]+', '_', cat))
  ggsave(file.path(out, paste0(name, '.pdf')), p, width = 10, height = 6); ggsave(file.path(out, paste0(name, '.png')), p, width = 10, height = 6, dpi = 180)
}
