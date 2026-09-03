#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(ggrepel)
})

set.seed(122383)

input_file <- file.path(
  "data", "public_reference", "GSE122383",
  "GSE122383_rnaseq_time_course_fpkm.xlsx"
)
out_dir <- file.path("output", "psc_kc_reference", "GSE122383")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- as.data.frame(read_excel(input_file))
raw <- raw[!is.na(raw$Gene) & nzchar(raw$Gene), ]
sample_cols <- c(
  "H9a", "H9b", "D7a", "D7b", "D14a", "D14b",
  "D21a", "D21b", "D45b", "D45c"
)
stopifnot(all(sample_cols %in% names(raw)))

fpkm <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- raw$Gene
if (anyDuplicated(rownames(fpkm))) {
  fpkm <- rowsum(fpkm, group = rownames(fpkm), reorder = FALSE)
}
log_expr <- log2(fpkm + 1)

stage <- c("D0", "D0", "D7", "D7", "D14", "D14", "D21", "D21", "D45", "D45")
day <- c(0, 0, 7, 7, 14, 14, 21, 21, 45, 45)
replicate <- rep(1:2, 5)
meta <- data.frame(
  sample = sample_cols,
  stage = stage,
  day = day,
  replicate = replicate,
  psc_derived = TRUE,
  include_in_calendar_validation = TRUE,
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

con <- gzfile(file.path(out_dir, "normalized_log2fpkm.tsv.gz"), "wt")
write.table(log_expr, con, sep = "\t", quote = FALSE, col.names = NA)
close(con)

stage_levels <- c("D0", "D7", "D14", "D21", "D45")
stage_days <- c(D0 = 0, D7 = 7, D14 = 14, D21 = 21, D45 = 45)
stage_mean <- sapply(stage_levels, function(s) {
  rowMeans(log_expr[, meta$sample[meta$stage == s], drop = FALSE])
})
write.table(
  stage_mean, file.path(out_dir, "stage_mean_log2fpkm.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

modules <- list(
  Pluripotency = c("POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"),
  Surface_ectoderm = c("KRT8", "KRT18", "KRT19", "TFAP2A", "TFAP2C", "EPCAM"),
  Epidermal_commitment = c("TP63", "KRT5", "KRT14", "ITGA6", "ITGB4", "COL17A1"),
  Early_maturation = c("KRT1", "KRT10", "IVL", "TGM1", "DSG1", "CLDN1"),
  Terminal_barrier = c("FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B"),
  Off_target = c("T", "SOX1", "PAX6", "VIM", "COL1A1")
)

gene_z <- t(scale(t(log_expr)))
gene_z[!is.finite(gene_z)] <- 0
module_scores <- do.call(rbind, lapply(names(modules), function(module) {
  genes <- intersect(modules[[module]], rownames(gene_z))
  data.frame(
    sample = colnames(gene_z),
    module = module,
    score = colMeans(gene_z[genes, , drop = FALSE]),
    genes_used = length(genes),
    stringsAsFactors = FALSE
  )
}))
module_scores <- merge(module_scores, meta, by = "sample")
write.csv(module_scores, file.path(out_dir, "module_scores.csv"), row.names = FALSE)

p_modules <- ggplot(
  subset(module_scores, module != "Off_target"),
  aes(day, score, color = module)
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2.7) +
  geom_point(alpha = 0.45, size = 1.6) +
  scale_x_continuous(breaks = unname(stage_days)) +
  labs(
    title = "GSE122383 PSC-derived epidermal time course",
    subtitle = "Gene-wise standardized log2(FPKM+1); two replicates per stage",
    x = "Differentiation day", y = "Module score", color = "Module"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "module_dynamics.png"), p_modules,
       width = 10, height = 6.5, dpi = 240)

variable <- apply(log_expr, 1, var)
top <- names(sort(variable, decreasing = TRUE))[seq_len(min(2000, length(variable)))]
pca <- prcomp(t(log_expr[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- cbind(meta, PC1 = pca$x[meta$sample, 1], PC2 = pca$x[meta$sample, 2])
write.csv(pca_df, file.path(out_dir, "pca_coordinates.csv"), row.names = FALSE)

centroid <- aggregate(cbind(PC1, PC2) ~ day, pca_df, mean)
centroid <- centroid[order(centroid$day), ]
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = factor(day))) +
  geom_path(
    data = centroid, aes(PC1, PC2, group = 1), inherit.aes = FALSE,
    color = "grey45", linewidth = 0.8,
    arrow = arrow(length = grid::unit(0.15, "cm"))
  ) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  labs(
    title = "GSE122383 PSC-derived trajectory",
    subtitle = "PCA of 2,000 most variable genes",
    x = sprintf("PC1 (%.1f%%)", pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_var[2]), color = "Day"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "pca_reference_trajectory.png"), p_pca,
       width = 8.5, height = 6.2, dpi = 240)

core_genes <- c(
  "POU5F1", "NANOG", "TFAP2A", "KRT18", "KRT19",
  "TP63", "KRT5", "KRT14", "ITGA6", "KRT1", "IVL",
  "SPRR1B", "ABCA12", "DSG1", "TGM1"
)
core_genes <- intersect(core_genes, rownames(log_expr))

build_template <- function(train_samples) {
  anchors <- sapply(stage_levels, function(s) {
    sample_name <- intersect(train_samples, meta$sample[meta$stage == s])
    stopifnot(length(sample_name) == 1)
    log_expr[core_genes, sample_name]
  })
  center <- rowMeans(anchors)
  spread <- apply(anchors, 1, sd)
  keep <- is.finite(spread) & spread > 0.15
  anchors <- anchors[keep, , drop = FALSE]
  center <- center[keep]
  spread <- spread[keep]
  daily <- t(apply(anchors, 1, function(x) {
    approx(unname(stage_days), x, xout = 0:45, rule = 2)$y
  }))
  daily <- sweep(sweep(daily, 1, center, "-"), 1, spread, "/")
  list(daily = daily, center = center, spread = spread, genes = rownames(anchors))
}

pred_rows <- list()
idx <- 1
for (train_rep in 1:2) {
  train_samples <- meta$sample[meta$replicate == train_rep]
  test_samples <- meta$sample[meta$replicate != train_rep]
  template <- build_template(train_samples)
  for (query in test_samples) {
    query_z <- (log_expr[template$genes, query] - template$center) / template$spread
    rmse <- sqrt(colMeans((template$daily - query_z)^2))
    pred_day <- which.min(rmse) - 1
    pred_rows[[idx]] <- data.frame(
      training_replicate = train_rep,
      validation_replicate = 3 - train_rep,
      sample = query,
      actual_stage = meta$stage[match(query, meta$sample)],
      actual_day = meta$day[match(query, meta$sample)],
      predicted_day = pred_day,
      absolute_error_days = abs(pred_day - meta$day[match(query, meta$sample)]),
      best_rmse = min(rmse),
      n_genes = length(template$genes),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}
pred <- do.call(rbind, pred_rows)
write.csv(pred, file.path(out_dir, "internal_cross_replicate_predictions.csv"),
          row.names = FALSE)

metrics <- data.frame(
  n_validation_samples = nrow(pred),
  mean_absolute_error_days = mean(pred$absolute_error_days),
  median_absolute_error_days = median(pred$absolute_error_days),
  within_2_days_fraction = mean(pred$absolute_error_days <= 2),
  within_4_days_fraction = mean(pred$absolute_error_days <= 4),
  maximum_absolute_error_days = max(pred$absolute_error_days),
  n_core_genes_available = length(core_genes),
  n_core_genes_used_min = min(pred$n_genes),
  n_core_genes_used_max = max(pred$n_genes)
)
write.csv(metrics, file.path(out_dir, "internal_cross_replicate_metrics.csv"),
          row.names = FALSE)

gse120 <- read.delim(
  file.path("output", "psc_kc_reference", "GSE120107", "stage_mean_log2cpm.tsv"),
  check.names = FALSE, row.names = 1
)
gse120_days <- c(D0 = 0, D4 = 4, D7 = 7, D15 = 15, D30 = 30)
common <- intersect(core_genes, rownames(gse120))
direction <- do.call(rbind, lapply(common, function(gene) {
  rho120 <- suppressWarnings(cor(
    as.numeric(gse120[gene, names(gse120_days)]), unname(gse120_days),
    method = "spearman"
  ))
  rho122 <- suppressWarnings(cor(
    as.numeric(stage_mean[gene, names(stage_days)]), unname(stage_days),
    method = "spearman"
  ))
  data.frame(
    gene = gene,
    spearman_GSE120107 = rho120,
    spearman_GSE122383 = rho122,
    direction_concordant = sign(rho120) == sign(rho122),
    stringsAsFactors = FALSE
  )
}))
write.csv(direction, file.path(out_dir, "cross_study_core_gene_direction.csv"),
          row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("GSE122383 analysis written to: ", out_dir)
