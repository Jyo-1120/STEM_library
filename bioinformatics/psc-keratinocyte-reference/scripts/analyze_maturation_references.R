#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

output_dir <- "output/psc_kc_reference/maturation_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

marker_modules <- list(
  Basal_identity = c("TP63", "KRT5", "KRT14", "COL17A1", "ITGB4", "ITGA6"),
  Early_maturation = c("KRT1", "KRT10", "IVL", "DSG1", "TGM1", "CLDN1"),
  Terminal_barrier = c("FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B")
)

marker_table <- do.call(
  rbind,
  lapply(names(marker_modules), function(module) {
    data.frame(module = module, gene = marker_modules[[module]], stringsAsFactors = FALSE)
  })
)

ensembl_to_symbol <- function(ensembl_ids) {
  clean_ids <- sub("\\..*$", "", ensembl_ids)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = unique(clean_ids),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(symbols[clean_ids])
}

collapse_counts_by_symbol <- function(count_matrix, ensembl_ids) {
  symbols <- ensembl_to_symbol(ensembl_ids)
  keep <- !is.na(symbols) & nzchar(symbols)
  collapsed <- rowsum(count_matrix[keep, , drop = FALSE], group = symbols[keep], reorder = FALSE)
  storage.mode(collapsed) <- "integer"
  collapsed
}

normalize_counts <- function(count_matrix, group) {
  y <- DGEList(counts = count_matrix, group = group)
  keep <- filterByExpr(y, group = group, min.count = 10)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  list(y = y, log_cpm = cpm(y, log = TRUE, prior.count = 2))
}

gene_zscore <- function(expression_matrix) {
  z <- t(scale(t(expression_matrix)))
  z[is.na(z)] <- 0
  z
}

make_module_scores <- function(expression_matrix, sample_meta, dataset_name) {
  result <- do.call(
    rbind,
    lapply(names(marker_modules), function(module) {
      genes <- intersect(marker_modules[[module]], rownames(expression_matrix))
      z <- gene_zscore(expression_matrix[genes, , drop = FALSE])
      data.frame(
        dataset = dataset_name,
        sample = colnames(expression_matrix),
        module = module,
        score = colMeans(z),
        genes_used = length(genes),
        stringsAsFactors = FALSE
      )
    })
  )
  merge(result, sample_meta, by = "sample", all.x = TRUE)
}

safe_gene_effect <- function(expression_matrix, gene, group_a, group_b, sample_groups) {
  if (!gene %in% rownames(expression_matrix)) return(NA_real_)
  mean(expression_matrix[gene, sample_groups == group_a]) -
    mean(expression_matrix[gene, sample_groups == group_b])
}

safe_gene_cor <- function(expression_matrix, gene, time_values) {
  if (!gene %in% rownames(expression_matrix)) return(NA_real_)
  suppressWarnings(cor(as.numeric(expression_matrix[gene, ]), time_values, method = "spearman"))
}

# -----------------------------------------------------------------------------
# GSE98483: primary keratinocyte D0/D2/D4/D7, two donor lines
# STAR ReadsPerGene columns are gene, unstranded, strand-1, strand-2. We use the
# unstranded column because the GEO record does not declare a counting strand.
# -----------------------------------------------------------------------------

gse984_files <- sort(list.files(
  "data/public_reference/GSE98483/extracted",
  pattern = "RNASeq.*ReadsPerGene.*gz$",
  full.names = TRUE
))
stopifnot(length(gse984_files) == 8)

gse984_list <- lapply(gse984_files, function(file_name) {
  tab <- read.delim(file_name, header = FALSE, skip = 4, stringsAsFactors = FALSE)
  data.frame(ensembl = tab[[1]], count = tab[[2]], stringsAsFactors = FALSE)
})
gse984_ensembl <- gse984_list[[1]]$ensembl
stopifnot(all(vapply(gse984_list, function(x) identical(x$ensembl, gse984_ensembl), logical(1))))

gse984_counts_ensembl <- do.call(cbind, lapply(gse984_list, `[[`, "count"))
gse984_sample_names <- sub("_ReadsPerGene.*$", "", basename(gse984_files))
colnames(gse984_counts_ensembl) <- gse984_sample_names
gse984_counts <- collapse_counts_by_symbol(gse984_counts_ensembl, gse984_ensembl)

gse984_meta <- data.frame(
  sample = gse984_sample_names,
  day = as.numeric(sub(".*-day([0-9]+).*", "\\1", gse984_sample_names)),
  donor = ifelse(grepl("PKC19", gse984_sample_names), "PKC19", "Dombi23"),
  stringsAsFactors = FALSE
)
gse984_meta$stage <- paste0("D", gse984_meta$day)
gse984_meta <- gse984_meta[order(gse984_meta$day, gse984_meta$donor), ]
gse984_counts <- gse984_counts[, gse984_meta$sample, drop = FALSE]
gse984_norm <- normalize_counts(gse984_counts, factor(gse984_meta$stage, levels = c("D0", "D2", "D4", "D7")))
gse984_logcpm <- gse984_norm$log_cpm

write.csv(gse984_meta, file.path(output_dir, "GSE98483_sample_metadata.csv"), row.names = FALSE)
gse984_connection <- gzfile(file.path(output_dir, "GSE98483_log2cpm.tsv.gz"), "wt")
write.table(gse984_logcpm, gse984_connection, sep = "\t", quote = FALSE, col.names = NA)
close(gse984_connection)

gse984_module_scores <- make_module_scores(gse984_logcpm, gse984_meta, "GSE98483")

gse984_variance <- apply(gse984_logcpm, 1, var)
gse984_top <- names(sort(gse984_variance, decreasing = TRUE))[seq_len(min(2000, length(gse984_variance)))]
gse984_pca <- prcomp(t(gse984_logcpm[gse984_top, , drop = FALSE]), center = TRUE)
gse984_pca_var <- 100 * gse984_pca$sdev^2 / sum(gse984_pca$sdev^2)
gse984_pca_df <- cbind(
  gse984_meta,
  PC1 = gse984_pca$x[gse984_meta$sample, 1],
  PC2 = gse984_pca$x[gse984_meta$sample, 2]
)
gse984_centroids <- aggregate(cbind(PC1, PC2) ~ day, data = gse984_pca_df, mean)
gse984_centroids <- gse984_centroids[order(gse984_centroids$day), ]

p_gse984_pca <- ggplot(gse984_pca_df, aes(PC1, PC2, color = factor(day), shape = donor)) +
  geom_path(
    data = gse984_centroids,
    aes(x = PC1, y = PC2, group = 1),
    inherit.aes = FALSE,
    color = "#555555", linewidth = 0.8,
    arrow = arrow(length = grid::unit(0.15, "cm"))
  ) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = paste0(donor, "-D", day)),
    size = 3.1, show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_color_brewer(palette = "YlOrRd", name = "Day") +
  labs(
    title = "GSE98483 primary-keratinocyte maturation",
    subtitle = "Growth-factor depletion/contact inhibition; two donor lines",
    x = sprintf("PC1 (%.1f%%)", gse984_pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", gse984_pca_var[2]),
    shape = "Donor"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "GSE98483_maturation_pca.png"), p_gse984_pca, width = 9, height = 6.5, dpi = 300)

# -----------------------------------------------------------------------------
# GSE73305: primary keratinocyte Ca2+ differentiation UD/D1-D5, no replicates
# -----------------------------------------------------------------------------

gse733_raw <- read.delim(
  "data/public_reference/GSE73305/GSE73305_All_genes_plus_gene_names.txt.gz",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gse733_expression <- as.matrix(gse733_raw[, c("UD", "D1", "D2", "D3", "D4", "D5")])
storage.mode(gse733_expression) <- "numeric"
gse733_symbols <- gse733_raw$Symbol
valid_733 <- !is.na(gse733_symbols) & nzchar(gse733_symbols)
gse733_expression <- rowsum(gse733_expression[valid_733, , drop = FALSE], gse733_symbols[valid_733], reorder = FALSE)
gse733_expression <- log2(gse733_expression + 1)
gse733_meta <- data.frame(
  sample = colnames(gse733_expression),
  day = 0:5,
  donor = "single_series",
  stage = colnames(gse733_expression),
  stringsAsFactors = FALSE
)
gse733_module_scores <- make_module_scores(gse733_expression, gse733_meta, "GSE73305")

write.csv(gse733_meta, file.path(output_dir, "GSE73305_sample_metadata.csv"), row.names = FALSE)
write.table(
  gse733_expression,
  file.path(output_dir, "GSE73305_log2_expression.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

# -----------------------------------------------------------------------------
# GSE287810: iPSC-derived KC in CnT30 versus KSFM
# -----------------------------------------------------------------------------

gse287_raw <- read.delim(
  "data/public_reference/GSE287810/GSE287810_LM_bulkRNAseq_iKC_rawcounts.tsv.gz",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gse287_counts_ensembl <- as.matrix(gse287_raw[, -1, drop = FALSE])
storage.mode(gse287_counts_ensembl) <- "integer"
gse287_counts <- collapse_counts_by_symbol(gse287_counts_ensembl, gse287_raw$gene)
gse287_meta <- data.frame(
  sample = colnames(gse287_counts),
  medium = ifelse(grepl("^cnt30", colnames(gse287_counts)), "CnT30", "KSFM"),
  batch = ifelse(grepl("_52_", colnames(gse287_counts)), "Diff52", "Diff6"),
  stringsAsFactors = FALSE
)
gse287_meta$day <- ifelse(gse287_meta$sample == "cnt30_52_1", 32, 29)
gse287_meta$stage <- gse287_meta$medium
gse287_norm <- normalize_counts(gse287_counts, factor(gse287_meta$medium))
gse287_logcpm <- gse287_norm$log_cpm
gse287_module_scores <- make_module_scores(gse287_logcpm, gse287_meta, "GSE287810")

write.csv(gse287_meta, file.path(output_dir, "GSE287810_sample_metadata.csv"), row.names = FALSE)
gse287_connection <- gzfile(file.path(output_dir, "GSE287810_log2cpm.tsv.gz"), "wt")
write.table(gse287_logcpm, gse287_connection, sep = "\t", quote = FALSE, col.names = NA)
close(gse287_connection)

design_287 <- model.matrix(~batch + medium, data = gse287_meta)
rownames(design_287) <- gse287_meta$sample
y287 <- estimateDisp(gse287_norm$y, design_287, robust = TRUE)
fit287 <- glmQLFit(y287, design_287, robust = TRUE)
coef_medium <- grep("medium", colnames(design_287))
qlf287 <- glmQLFTest(fit287, coef = coef_medium)
de287 <- topTags(qlf287, n = Inf, sort.by = "PValue")$table
de287$gene <- rownames(de287)
write.csv(de287[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")], file.path(output_dir, "GSE287810_KSFM_vs_CnT30_DE.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Module dynamics in two primary-KC maturation references
# -----------------------------------------------------------------------------

primary_module_scores <- rbind(gse984_module_scores, gse733_module_scores)
primary_module_scores$day <- as.numeric(primary_module_scores$day)
write.csv(primary_module_scores, file.path(output_dir, "primary_KC_maturation_module_scores.csv"), row.names = FALSE)

module_colors <- c(
  Basal_identity = "#E5A11A",
  Early_maturation = "#F58518",
  Terminal_barrier = "#B279A2"
)

p_primary_modules <- ggplot(primary_module_scores, aes(day, score, color = module, group = interaction(module, donor))) +
  geom_line(alpha = 0.65, linewidth = 0.8) +
  geom_point(size = 2.5) +
  stat_summary(aes(group = module), fun = mean, geom = "line", linewidth = 1.2) +
  scale_color_manual(values = module_colors) +
  facet_wrap(~dataset, scales = "free_x") +
  labs(
    title = "Primary-keratinocyte maturation module dynamics",
    subtitle = "Within-study gene-wise standardized scores; studies are not merged at count level",
    x = "Days after differentiation induction", y = "Module score", color = "Module"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "primary_KC_maturation_modules.png"), p_primary_modules, width = 11, height = 6.5, dpi = 300)

p_287_modules <- ggplot(gse287_module_scores, aes(medium, score, color = module)) +
  geom_point(position = position_jitter(width = 0.08), size = 2.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.7) +
  facet_wrap(~module, nrow = 1) +
  scale_color_manual(values = module_colors) +
  labs(
    title = "GSE287810 iPSC-KC module comparison",
    subtitle = "CnT30 versus KSFM at approximately D29-D32",
    x = "Keratinocyte medium", y = "Within-dataset module score"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "GSE287810_medium_module_comparison.png"), p_287_modules, width = 10.5, height = 5.2, dpi = 300)

# -----------------------------------------------------------------------------
# Cross-dataset marker validation
# -----------------------------------------------------------------------------

gse120_means <- as.matrix(read.delim(
  "output/psc_kc_reference/GSE120107/stage_mean_log2cpm.tsv",
  row.names = 1,
  check.names = FALSE
))

validation_rows <- lapply(seq_len(nrow(marker_table)), function(i) {
  gene <- marker_table$gene[i]
  module <- marker_table$module[i]

  effect_120_d30_d15 <- if (gene %in% rownames(gse120_means)) gse120_means[gene, "D30"] - gse120_means[gene, "D15"] else NA_real_
  effect_120_primary_d30 <- if (gene %in% rownames(gse120_means)) gse120_means[gene, "Primary_KC"] - gse120_means[gene, "D30"] else NA_real_
  effect_984 <- safe_gene_effect(gse984_logcpm, gene, "D7", "D0", gse984_meta$stage)
  rho_984 <- safe_gene_cor(gse984_logcpm, gene, gse984_meta$day)
  effect_733 <- safe_gene_effect(gse733_expression, gene, "D5", "UD", gse733_meta$stage)
  rho_733 <- safe_gene_cor(gse733_expression, gene, gse733_meta$day)
  # edgeR coefficient is KSFM minus CnT30; reverse the sign here so all output
  # columns consistently report CnT30 minus KSFM while controlling for batch.
  effect_287 <- if (gene %in% de287$gene) -de287$logFC[match(gene, de287$gene)] else NA_real_

  maturation_support <- sum(c(effect_120_primary_d30, effect_984, effect_733) >= 0.5, na.rm = TRUE)
  basal_support <- sum(c(effect_120_d30_d15, effect_287) >= 0.5, na.rm = TRUE)

  recommendation <- if (module %in% c("Early_maturation", "Terminal_barrier")) {
    if (maturation_support == 3) "direction-consistent maturation candidate"
    else if (maturation_support == 2) "partial support / context-dependent"
    else "low or inconsistent maturation support"
  } else {
    if (basal_support == 2) "direction-consistent basal iKC candidate"
    else if (basal_support == 1) "partial basal support"
    else "low or inconsistent basal support"
  }

  priority_for_qpcr <- if (module %in% c("Early_maturation", "Terminal_barrier")) {
    maturation_support == 3 && !is.na(rho_984) && rho_984 >= 0.5 && !is.na(rho_733) && rho_733 >= 0.5
  } else {
    basal_support == 2 && !is.na(effect_984) && effect_984 <= -0.3 && !is.na(effect_733) && effect_733 <= -0.3
  }

  data.frame(
    module = module,
    gene = gene,
    GSE120107_D30_vs_D15 = effect_120_d30_d15,
    GSE120107_primaryKC_vs_D30 = effect_120_primary_d30,
    GSE98483_D7_vs_D0 = effect_984,
    GSE98483_spearman_day = rho_984,
    GSE73305_D5_vs_UD = effect_733,
    GSE73305_spearman_day = rho_733,
    GSE287810_CnT30_vs_KSFM = effect_287,
    maturation_support_0_to_3 = maturation_support,
    basal_support_0_to_2 = basal_support,
    priority_for_qpcr_bulk_validation = priority_for_qpcr,
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
})

validation <- do.call(rbind, validation_rows)
validation <- validation[order(validation$module, -validation$maturation_support_0_to_3, -validation$basal_support_0_to_2), ]
write.csv(validation, file.path(output_dir, "cross_dataset_marker_validation.csv"), row.names = FALSE)

effect_matrix <- as.matrix(validation[, c(
  "GSE120107_D30_vs_D15",
  "GSE120107_primaryKC_vs_D30",
  "GSE98483_D7_vs_D0",
  "GSE73305_D5_vs_UD",
  "GSE287810_CnT30_vs_KSFM"
)])
rownames(effect_matrix) <- validation$gene
effect_matrix[effect_matrix > 6] <- 6
effect_matrix[effect_matrix < -6] <- -6
effect_colnames <- c("PSC-KC D30-D15", "PrimaryKC-D30", "PrimaryKC D7-D0", "Ca D5-UD", "CnT30-KSFM")
colnames(effect_matrix) <- effect_colnames
effect_annotation <- data.frame(Module = validation$module)
rownames(effect_annotation) <- validation$gene

png(file.path(output_dir, "cross_dataset_marker_effects_heatmap.png"), width = 2500, height = 2600, res = 300)
pheatmap(
  effect_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = effect_annotation,
  annotation_colors = list(Module = module_colors),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(121),
  breaks = seq(-6, 6, length.out = 122),
  na_col = "#D9D9D9",
  main = "Cross-dataset marker effects (log2-scale differences, capped at +/-6)",
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 8
)
dev.off()

consistent_maturation <- validation$gene[validation$recommendation == "direction-consistent maturation candidate"]
priority_maturation <- validation$gene[
  validation$priority_for_qpcr_bulk_validation & validation$module %in% c("Early_maturation", "Terminal_barrier")
]
priority_basal <- validation$gene[
  validation$priority_for_qpcr_bulk_validation & validation$module == "Basal_identity"
]

summary_lines <- c(
  "# 후기 keratinocyte maturation reference: 교차 데이터 검증",
  "",
  "## 분석 범위",
  "",
  "- GSE98483: primary KC D0/D2/D4/D7, 2 donor lines",
  "- GSE73305: calcium-induced primary KC UD/D1-D5, biological replicate 없음",
  "- GSE287810: iPSC-KC CnT30 대 KSFM, 약 D29-D32",
  "- GSE120107: 이전 분석의 D15/D30/primary KC stage mean",
  "",
  "각 연구는 별도로 정규화했으며 count matrix를 직접 병합하지 않았다.",
  "",
  "## qPCR marker 교차 검증 결과",
  "",
  sprintf("- 세 maturation 비교에서 증가 방향이 일치한 후보: %s", if (length(consistent_maturation)) paste(consistent_maturation, collapse = ", ") else "없음"),
  sprintf("- 두 primary KC 자료에서도 시간에 따라 단조롭게 증가한 우선 qPCR 후보: %s", if (length(priority_maturation)) paste(priority_maturation, collapse = ", ") else "없음"),
  sprintf("- iKC에서 증가하고 primary KC maturation에서 감소한 basal 축 우선 후보: %s", if (length(priority_basal)) paste(priority_basal, collapse = ", ") else "없음"),
  "",
  "## 해석 원칙",
  "",
  "- Basal marker와 maturation marker를 함께 측정해야 증식성 iKC와 성숙 KC를 구분할 수 있다.",
  "- GSE73305는 반복이 없으므로 통계 검증이 아니라 방향성 확인에만 사용한다.",
  "- CnT30-KSFM 비교는 성숙 시간축이 아니라 iKC 배양 조건에 따른 keratinocyte identity 보강 자료이다.",
  "- 최종 qPCR panel은 이후 scRNA-seq에서 epithelial population 특이성을 확인한 뒤 확정한다."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary_ko.md"))

message("Maturation-reference analysis completed: ", output_dir)
