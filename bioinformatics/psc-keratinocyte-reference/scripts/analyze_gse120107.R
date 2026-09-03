#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

input_file <- "data/public_reference/GSE120107/GSE120107_normal_RNA-seq_counts.txt.gz"
output_dir <- "output/psc_kc_reference/GSE120107"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stage_levels <- c("D0", "D4", "D7", "D15", "D30", "Primary_KC")
stage_order <- setNames(seq_along(stage_levels) - 1, stage_levels)
stage_colors <- c(
  D0 = "#4C78A8", D4 = "#72B7B2", D7 = "#54A24B",
  D15 = "#F2CF5B", D30 = "#F58518", Primary_KC = "#B279A2"
)

counts_df <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(counts_df)[1] == "Genes")
stopifnot(!anyDuplicated(counts_df$Genes))

counts <- as.matrix(counts_df[, -1, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- counts_df$Genes

infer_stage <- function(sample_name) {
  if (grepl("_d0_", sample_name)) return("D0")
  if (grepl("_d4_", sample_name)) return("D4")
  if (grepl("_d7_", sample_name)) return("D7")
  if (grepl("_d15_", sample_name)) return("D15")
  if (grepl("_d30_", sample_name)) return("D30")
  if (grepl("_KRT_", sample_name)) return("Primary_KC")
  stop("Could not infer stage from sample name: ", sample_name)
}

sample_meta <- data.frame(
  sample = colnames(counts),
  stage = vapply(colnames(counts), infer_stage, character(1)),
  stringsAsFactors = FALSE
)
sample_meta$stage <- factor(sample_meta$stage, levels = stage_levels)
sample_meta$stage_order <- unname(stage_order[as.character(sample_meta$stage)])
sample_meta$replicate <- ave(
  seq_len(nrow(sample_meta)), sample_meta$stage,
  FUN = function(x) seq_along(x)
)
rownames(sample_meta) <- sample_meta$sample

sample_meta <- sample_meta[order(sample_meta$stage, sample_meta$replicate), ]
counts <- counts[, sample_meta$sample, drop = FALSE]

y_all <- DGEList(counts = counts, group = sample_meta$stage)
sample_meta$library_size <- colSums(counts)
sample_meta$detected_genes_raw_count_gt_0 <- colSums(counts > 0)

keep <- filterByExpr(y_all, group = sample_meta$stage, min.count = 10)
y <- y_all[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y, method = "TMM")
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
sample_meta$tmm_norm_factor <- y$samples$norm.factors
sample_meta$effective_library_size <- y$samples$lib.size * y$samples$norm.factors
sample_meta$detected_genes_cpm_ge_1 <- colSums(cpm(y) >= 1)

write.csv(sample_meta, file.path(output_dir, "sample_qc.csv"), row.names = FALSE)
normalized_connection <- gzfile(file.path(output_dir, "normalized_log2cpm.tsv.gz"), "wt")
write.table(log_cpm, normalized_connection, sep = "\t", quote = FALSE, col.names = NA)
close(normalized_connection)

gene_variance <- apply(log_cpm, 1, var)
top_variable <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(2000, length(gene_variance)))]
pca <- prcomp(t(log_cpm[top_variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- cbind(
  sample_meta,
  PC1 = pca$x[sample_meta$sample, 1],
  PC2 = pca$x[sample_meta$sample, 2]
)
write.csv(pca_df, file.path(output_dir, "pca_coordinates.csv"), row.names = FALSE)

centroids <- aggregate(cbind(PC1, PC2) ~ stage + stage_order, data = pca_df, FUN = mean)
centroids <- centroids[order(centroids$stage_order), ]

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = stage)) +
  geom_path(
    data = centroids, aes(x = PC1, y = PC2, group = 1), inherit.aes = FALSE,
    color = "#555555", linewidth = 0.8, arrow = arrow(length = unit(0.15, "cm"))
  ) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = paste0(stage, "-R", replicate)),
    size = 3.2, show.legend = FALSE, max.overlaps = Inf,
    box.padding = 0.45, point.padding = 0.25, min.segment.length = 0
  ) +
  scale_color_manual(values = stage_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(
    title = "GSE120107 PSC-to-keratinocyte reference trajectory",
    subtitle = "TMM-normalized log2 CPM; PCA of 2,000 most variable genes",
    x = sprintf("PC1 (%.1f%%)", pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_var[2]),
    color = "Reference stage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(12, 24, 12, 12)
  )

ggsave(file.path(output_dir, "pca_reference_trajectory.png"), p_pca, width = 9, height = 6.5, dpi = 300)

sample_cor <- cor(log_cpm, method = "pearson")
annotation_col <- data.frame(Stage = sample_meta$stage)
rownames(annotation_col) <- sample_meta$sample
annotation_colors <- list(Stage = stage_colors)

png(file.path(output_dir, "sample_correlation_heatmap.png"), width = 2400, height = 2100, res = 300)
pheatmap(
  sample_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  main = "GSE120107 sample correlation",
  border_color = NA,
  fontsize = 8
)
dev.off()

marker_modules <- list(
  Pluripotency = c("POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"),
  Surface_ectoderm = c("KRT8", "KRT18", "KRT19", "TFAP2A", "TFAP2C", "EPCAM"),
  Epidermal_commitment = c("TP63", "KRT5", "KRT14", "ITGA6", "ITGB4", "COL17A1"),
  Early_maturation = c("KRT1", "KRT10", "IVL", "TGM1", "DSG1", "CLDN1"),
  Terminal_barrier = c("FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B"),
  Off_target = c("T", "SOX1", "PAX6", "VIM", "COL1A1")
)

marker_table <- do.call(
  rbind,
  lapply(names(marker_modules), function(module) {
    data.frame(module = module, gene = marker_modules[[module]], stringsAsFactors = FALSE)
  })
)
marker_table$detected <- marker_table$gene %in% rownames(log_cpm)
write.csv(marker_table, file.path(output_dir, "marker_detection.csv"), row.names = FALSE)

present_markers <- marker_table[marker_table$detected, , drop = FALSE]
marker_expression <- log_cpm[present_markers$gene, , drop = FALSE]
marker_z <- t(scale(t(marker_expression)))
marker_z[is.na(marker_z)] <- 0

annotation_row <- data.frame(Module = present_markers$module)
rownames(annotation_row) <- present_markers$gene
module_colors <- c(
  Pluripotency = "#4C78A8", Surface_ectoderm = "#72B7B2",
  Epidermal_commitment = "#ECA82C", Early_maturation = "#F58518",
  Terminal_barrier = "#B279A2", Off_target = "#9D755D"
)

png(file.path(output_dir, "stage_marker_heatmap.png"), width = 2600, height = 3000, res = 300)
pheatmap(
  marker_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = list(Module = module_colors, Stage = stage_colors),
  color = colorRampPalette(c("#2C3E75", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  main = "Stage-marker expression across PSC-to-keratinocyte differentiation",
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 7
)
dev.off()

write.table(
  marker_expression,
  file.path(output_dir, "marker_log2cpm.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

stage_means <- sapply(stage_levels, function(stage_name) {
  rowMeans(log_cpm[, sample_meta$sample[sample_meta$stage == stage_name], drop = FALSE])
})
write.table(
  stage_means,
  file.path(output_dir, "stage_mean_log2cpm.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

module_scores <- do.call(
  rbind,
  lapply(names(marker_modules), function(module) {
    genes <- intersect(marker_modules[[module]], rownames(log_cpm))
    gene_z <- t(scale(t(log_cpm[genes, , drop = FALSE])))
    gene_z[is.na(gene_z)] <- 0
    data.frame(
      sample = colnames(log_cpm),
      module = module,
      score = colMeans(gene_z),
      genes_used = length(genes),
      stringsAsFactors = FALSE
    )
  })
)
module_scores <- merge(module_scores, sample_meta[, c("sample", "stage", "stage_order", "replicate")], by = "sample")
module_scores$stage <- factor(module_scores$stage, levels = stage_levels)
write.csv(module_scores, file.path(output_dir, "module_scores.csv"), row.names = FALSE)

p_modules <- ggplot(module_scores[module_scores$module != "Off_target", ], aes(stage_order, score, color = module)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2.8) +
  geom_point(alpha = 0.45, size = 1.7) +
  scale_x_continuous(breaks = seq_along(stage_levels) - 1, labels = stage_levels) +
  scale_color_manual(values = module_colors) +
  labs(
    title = "Reference gene-module dynamics",
    subtitle = "Scores are means of gene-wise standardized log2 CPM values",
    x = "Reference stage", y = "Module score", color = "Module"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "module_dynamics.png"), p_modules, width = 9.5, height = 6.5, dpi = 300)

candidate_rows <- lapply(seq_len(nrow(present_markers)), function(i) {
  gene <- present_markers$gene[i]
  values <- log_cpm[gene, ]
  means <- stage_means[gene, ]
  data.frame(
    module = present_markers$module[i],
    gene = gene,
    peak_stage = names(which.max(means)),
    dynamic_range_log2cpm = max(means) - min(means),
    log2fc_D30_vs_D0 = means["D30"] - means["D0"],
    log2fc_primaryKC_vs_D0 = means["Primary_KC"] - means["D0"],
    spearman_with_stage = suppressWarnings(cor(values, sample_meta$stage_order, method = "spearman")),
    stringsAsFactors = FALSE
  )
})
candidate_summary <- do.call(rbind, candidate_rows)
candidate_summary <- candidate_summary[order(candidate_summary$module, -candidate_summary$dynamic_range_log2cpm), ]
write.csv(candidate_summary, file.path(output_dir, "qpcr_candidate_summary.csv"), row.names = FALSE)

housekeeping_candidates <- c("RPLP0", "HPRT1", "TBP", "PPIA", "GAPDH", "ACTB", "B2M")
housekeeping_candidates <- intersect(housekeeping_candidates, rownames(log_cpm))
housekeeping <- data.frame(
  gene = housekeeping_candidates,
  mean_log2cpm = rowMeans(log_cpm[housekeeping_candidates, , drop = FALSE]),
  sd_log2cpm = apply(log_cpm[housekeeping_candidates, , drop = FALSE], 1, sd),
  range_log2cpm = apply(log_cpm[housekeeping_candidates, , drop = FALSE], 1, function(x) diff(range(x))),
  stringsAsFactors = FALSE
)
housekeeping <- housekeeping[order(housekeeping$sd_log2cpm, housekeeping$range_log2cpm), ]
write.csv(housekeeping, file.path(output_dir, "housekeeping_stability.csv"), row.names = FALSE)

design <- model.matrix(~0 + stage, data = sample_meta)
colnames(design) <- stage_levels
y_de <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y_de, design, robust = TRUE)
contrast_strings <- c(
  D4_vs_D0 = "D4-D0",
  D7_vs_D4 = "D7-D4",
  D15_vs_D7 = "D15-D7",
  D30_vs_D15 = "D30-D15",
  PrimaryKC_vs_D30 = "Primary_KC-D30"
)

de_results <- do.call(
  rbind,
  lapply(names(contrast_strings), function(contrast_name) {
    contrast <- makeContrasts(contrasts = contrast_strings[[contrast_name]], levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    tab$gene <- rownames(tab)
    tab$contrast <- contrast_name
    tab[, c("contrast", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  })
)
de_connection <- gzfile(file.path(output_dir, "adjacent_stage_differential_expression.tsv.gz"), "wt")
write.table(de_results, de_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(de_connection)

de_summary <- do.call(
  rbind,
  lapply(names(contrast_strings), function(contrast_name) {
    tab <- de_results[de_results$contrast == contrast_name, ]
    data.frame(
      contrast = contrast_name,
      significant_up = sum(tab$FDR < 0.05 & tab$logFC >= 1),
      significant_down = sum(tab$FDR < 0.05 & tab$logFC <= -1),
      tested_genes = nrow(tab),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(de_summary, file.path(output_dir, "transition_de_summary.csv"), row.names = FALSE)

top_transition_genes <- do.call(
  rbind,
  lapply(names(contrast_strings), function(contrast_name) {
    tab <- de_results[de_results$contrast == contrast_name & de_results$FDR < 0.05 & abs(de_results$logFC) >= 1, ]
    if (!nrow(tab)) return(NULL)
    do.call(rbind, lapply(c("up", "down"), function(direction) {
      direction_tab <- if (direction == "up") tab[tab$logFC > 0, ] else tab[tab$logFC < 0, ]
      if (!nrow(direction_tab)) return(NULL)
      direction_tab <- direction_tab[order(direction_tab$FDR, -abs(direction_tab$logFC)), ]
      direction_tab <- head(direction_tab, 10)
      direction_tab$direction <- direction
      direction_tab[, c("contrast", "direction", "gene", "logFC", "logCPM", "FDR")]
    }))
  })
)
write.csv(top_transition_genes, file.path(output_dir, "top_transition_genes.csv"), row.names = FALSE)

within_stage_cor <- unlist(lapply(stage_levels, function(stage_name) {
  samples <- sample_meta$sample[sample_meta$stage == stage_name]
  if (length(samples) == 2) cor(log_cpm[, samples[1]], log_cpm[, samples[2]]) else NA_real_
}))
names(within_stage_cor) <- stage_levels

pc_stage_cor <- c(
  PC1 = cor(pca_df$PC1, pca_df$stage_order, method = "spearman"),
  PC2 = cor(pca_df$PC2, pca_df$stage_order, method = "spearman")
)

missing_markers <- marker_table$gene[!marker_table$detected]
summary_lines <- c(
  "# GSE120107 initial reference-map summary",
  "",
  sprintf("- Input samples: %d", ncol(counts)),
  sprintf("- Input genes: %s", format(nrow(counts), big.mark = ",")),
  sprintf("- Genes retained after expression filtering: %s", format(nrow(log_cpm), big.mark = ",")),
  sprintf("- Library-size range: %s to %s reads", format(min(sample_meta$library_size), big.mark = ","), format(max(sample_meta$library_size), big.mark = ",")),
  sprintf("- PC1/PC2 explained variance: %.1f%% / %.1f%%", pca_var[1], pca_var[2]),
  sprintf("- Spearman correlation with ordered reference stage: PC1 %.3f, PC2 %.3f", pc_stage_cor[1], pc_stage_cor[2]),
  sprintf("- Mean within-stage replicate correlation: %.4f", mean(within_stage_cor, na.rm = TRUE)),
  sprintf("- Within-stage replicate correlations: %s", paste(sprintf("%s=%.3f", names(within_stage_cor), within_stage_cor), collapse = ", ")),
  sprintf("- Most stable tested housekeeping genes: %s", paste(head(housekeeping$gene, 3), collapse = ", ")),
  sprintf("- Missing planned markers after filtering: %s", if (length(missing_markers)) paste(missing_markers, collapse = ", ") else "none"),
  sprintf(
    "- Differential genes by transition (up/down; FDR < 0.05 and |log2FC| >= 1): %s",
    paste(sprintf("%s=%d/%d", de_summary$contrast, de_summary$significant_up, de_summary$significant_down), collapse = ", ")
  ),
  "",
  "## Initial biological interpretation",
  "",
  "- The pluripotency module falls sharply from D0 and is largely suppressed by D15.",
  "- The surface-ectoderm/simple-epithelium module is highest at D4-D7 and decreases after epidermal commitment.",
  "- TP63/KRT5/KRT14/COL17A1-associated epidermal commitment rises strongly from D15 to D30.",
  "- Terminal/barrier markers remain substantially lower at D30 than in primary keratinocytes, indicating that D30 is closer to an immature or proliferative basal-KC state than to a fully mature epidermal state.",
  "- RPLP0, PPIA, and TBP are provisional qPCR reference-gene candidates in this dataset; stability must be retested in the lab's own samples.",
  "",
  "The PCA trajectory is exploratory. Public differentiation day and biological state are not interchangeable; future lab qPCR samples require within-lab calibration."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.md"))

message("Analysis completed: ", output_dir)
