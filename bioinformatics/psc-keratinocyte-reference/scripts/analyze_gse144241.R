#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
  library(ggplot2)
})

input_dir <- file.path("data", "public_reference", "GSE144241")
out_dir <- file.path("output", "psc_kc_reference", "GSE144241")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gset <- getGEO(
  filename = file.path(input_dir, "GSE144241_series_matrix.txt.gz"),
  getGPL = FALSE, AnnotGPL = FALSE
)
probe_expr <- exprs(gset)
colnames(probe_expr) <- c("D0", "D1", "D4", "D6", "D8", "D11", "D26")

annot_lines <- readLines(gzfile(file.path(input_dir, "GPL10558.annot.gz")))
header_line <- grep("^ID\\tGene title\\tGene symbol", annot_lines)[1]
footer_line <- grep("^!platform_table_end", annot_lines)[1]
annot_text <- annot_lines[header_line:(footer_line - 1)]
annot <- read.delim(
  textConnection(annot_text), check.names = FALSE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)
annot <- annot[, c("ID", "Gene symbol")]
names(annot)[2] <- "gene"
annot$gene <- sub(" ///.*$", "", annot$gene)
annot <- annot[nzchar(annot$gene), ]

common_probe <- intersect(rownames(probe_expr), annot$ID)
d <- data.frame(
  probe = common_probe,
  gene = annot$gene[match(common_probe, annot$ID)],
  probe_expr[common_probe, , drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
d$probe_variance <- apply(d[, colnames(probe_expr), drop = FALSE], 1, var)
d <- d[order(d$gene, -d$probe_variance), ]
d <- d[!duplicated(d$gene), ]
gene_expr <- as.matrix(d[, colnames(probe_expr), drop = FALSE])
rownames(gene_expr) <- d$gene

con <- gzfile(file.path(out_dir, "gene_level_quantile_normalized_expression.tsv.gz"), "wt")
write.table(gene_expr, con, sep = "\t", quote = FALSE, col.names = NA)
close(con)

meta <- data.frame(
  sample = colnames(gene_expr),
  day = c(0, 1, 4, 6, 8, 11, 26),
  replicate = 1,
  psc_derived = TRUE,
  training_eligible = FALSE,
  exclusion_reason = "one sample per time point",
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

modules <- list(
  Pluripotency = c("POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"),
  Surface_ectoderm = c("KRT8", "KRT18", "KRT19", "TFAP2A", "TFAP2C", "EPCAM"),
  Epidermal_commitment = c("TP63", "KRT5", "KRT14", "ITGA6", "ITGB4", "COL17A1"),
  Early_maturation = c("KRT1", "KRT10", "IVL", "TGM1", "DSG1", "CLDN1"),
  Terminal_barrier = c("FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B"),
  Off_target = c("T", "SOX1", "PAX6", "VIM", "COL1A1")
)

gene_z <- t(scale(t(gene_expr)))
gene_z[!is.finite(gene_z)] <- 0
module_scores <- do.call(rbind, lapply(names(modules), function(module) {
  genes <- intersect(modules[[module]], rownames(gene_z))
  data.frame(
    day = meta$day,
    sample = meta$sample,
    module = module,
    score = colMeans(gene_z[genes, , drop = FALSE]),
    genes_used = length(genes),
    stringsAsFactors = FALSE
  )
}))
write.csv(module_scores, file.path(out_dir, "module_scores.csv"), row.names = FALSE)

p <- ggplot(subset(module_scores, module != "Off_target"),
            aes(day, score, color = module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = meta$day) +
  labs(
    title = "GSE144241 hESC-to-keratinocyte time course",
    subtitle = "Quantile-normalized microarray; one sample per day",
    x = "Differentiation day", y = "Module score", color = "Module"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "module_dynamics.png"), p,
       width = 10, height = 6.5, dpi = 240)

core_genes <- c(
  "POU5F1", "NANOG", "TFAP2A", "KRT18", "KRT19",
  "TP63", "KRT5", "KRT14", "ITGA6", "KRT1", "IVL",
  "SPRR1B", "ABCA12", "DSG1", "TGM1"
)
present <- intersect(core_genes, rownames(gene_expr))
trend <- do.call(rbind, lapply(present, function(gene) {
  data.frame(
    gene = gene,
    spearman_with_day = suppressWarnings(cor(
      as.numeric(gene_expr[gene, ]), meta$day, method = "spearman"
    )),
    peak_day = meta$day[which.max(gene_expr[gene, ])],
    dynamic_range = diff(range(gene_expr[gene, ])),
    stringsAsFactors = FALSE
  )
}))
write.csv(trend, file.path(out_dir, "core_gene_trends.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("GSE144241 analysis written to: ", out_dir)
