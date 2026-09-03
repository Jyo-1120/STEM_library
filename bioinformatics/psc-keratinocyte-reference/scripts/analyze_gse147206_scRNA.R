#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(147206)

project_dir <- normalizePath(".")
data_dir <- file.path(project_dir, "data", "public_reference", "GSE147206", "extracted")
out_dir <- file.path(project_dir, "output", "psc_kc_reference", "GSE147206_scRNA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_info <- data.frame(
  sample = c("d6_WA25", "d6_DSP", "d29_WA25", "d29_DSP", "d48_WA25"),
  day = c("D6", "D6", "D29", "D29", "D48"),
  cell_line = c("WA25", "DSP", "WA25", "DSP", "WA25"),
  matrix_dir = file.path(
    data_dir,
    c("d6_WA25/filtered_feature_bc_matrix",
      "d6_DSP/filtered_feature_bc_matrix",
      "d29_WA25/filtered_gene_bc_matrices/GRCh38",
      "d29_DSP/filtered_gene_bc_matrices/GRCh38",
      "d48_WA25/filtered_gene_bc_matrices/GRCh38")
  ),
  stringsAsFactors = FALSE
)

read_sample <- function(i) {
  info <- sample_info[i, ]
  message("Reading ", info$sample)
  counts <- Read10X(info$matrix_dir, gene.column = 2)
  if (is.list(counts)) counts <- counts[[1]]
  colnames(counts) <- paste(info$sample, colnames(counts), sep = "_")
  obj <- CreateSeuratObject(
    counts = counts,
    project = info$sample,
    min.cells = 5,
    min.features = 200
  )
  obj$sample <- info$sample
  obj$day <- info$day
  obj$cell_line <- info$cell_line
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

objects_raw <- lapply(seq_len(nrow(sample_info)), read_sample)
names(objects_raw) <- sample_info$sample

qc_rows <- lapply(names(objects_raw), function(nm) {
  x <- objects_raw[[nm]]
  keep <- x$nFeature_RNA > 500 & x$nFeature_RNA < 7500 & x$percent.mt < 20
  data.frame(
    sample = nm,
    day = unique(x$day),
    cell_line = unique(x$cell_line),
    cells_after_minimum_filter = ncol(x),
    cells_after_qc = sum(keep),
    retained_fraction = mean(keep),
    median_features_before_qc = median(x$nFeature_RNA),
    median_features_after_qc = median(x$nFeature_RNA[keep]),
    median_mt_pct_after_qc = median(x$percent.mt[keep])
  )
})
qc_table <- do.call(rbind, qc_rows)
write.csv(qc_table, file.path(out_dir, "qc_cell_counts.csv"), row.names = FALSE)

objects <- lapply(objects_raw, function(x) {
  subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
})
rm(objects_raw)
invisible(gc())

integrate_day <- function(day_name, object_names, nfeatures, dims_use, resolution) {
  message("Integrating ", day_name)
  xs <- lapply(objects[object_names], function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    x
  })
  if (length(xs) > 1) {
    anchors <- FindIntegrationAnchors(object.list = xs, dims = dims_use, verbose = FALSE)
    out <- IntegrateData(anchorset = anchors, dims = dims_use, verbose = FALSE)
    DefaultAssay(out) <- "integrated"
  } else {
    out <- xs[[1]]
    DefaultAssay(out) <- "RNA"
  }
  out <- ScaleData(out, verbose = FALSE)
  out <- RunPCA(out, npcs = max(dims_use), verbose = FALSE)
  out <- FindNeighbors(out, dims = dims_use, verbose = FALSE)
  out <- FindClusters(out, resolution = resolution, verbose = FALSE)
  out <- RunUMAP(out, dims = dims_use, seed.use = 147206, verbose = FALSE)
  DefaultAssay(out) <- "RNA"
  out <- JoinLayers(out, assay = "RNA")
  out
}

day_objects <- list(
  D6 = integrate_day("D6", c("d6_WA25", "d6_DSP"), 2500, 1:40, 1.0),
  D29 = integrate_day("D29", c("d29_WA25", "d29_DSP"), 2000, 1:36, 0.8),
  D48 = integrate_day("D48", "d48_WA25", 2000, 1:30, 0.6)
)
rm(objects)
invisible(gc())

broad_modules <- list(
  Epithelial = c("EPCAM", "CDH1", "TACSTD2", "KRT8", "KRT18", "KRT19"),
  Mesenchymal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"),
  Neuroglial = c("SOX10", "PLP1", "S100B", "PAX6", "HES5", "NEUROD1"),
  Melanocyte = c("MLANA", "PMEL", "TYR", "DCT", "MITF"),
  Myocyte = c("MYOG", "MYOD1", "TNNT1", "DES")
)

epithelial_modules <- list(
  Basal = c("TP63", "KRT5", "KRT14", "COL17A1", "ITGA6", "ITGB4"),
  Periderm = c("KRT4", "KRT8", "KRT18", "KRT19", "KRT13", "SFN"),
  Intermediate = c("KRT1", "KRT10", "IVL", "DSG1", "TGM1", "CLDN1"),
  Terminal = c("FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B"),
  Cycling = c("MKI67", "TOP2A", "TYMS", "UBE2C", "CENPF")
)

add_scores <- function(obj, modules, prefix) {
  for (nm in names(modules)) {
    genes <- intersect(modules[[nm]], rownames(obj))
    obj <- AddModuleScore(obj, features = list(genes), name = paste0(prefix, nm), seed = 147206)
  }
  obj
}

annotate_clusters <- function(obj, day_name) {
  obj <- add_scores(obj, broad_modules, "broad_")
  obj <- add_scores(obj, epithelial_modules, "epi_")
  broad_cols <- paste0("broad_", names(broad_modules), "1")
  epi_cols <- paste0("epi_", names(epithelial_modules), "1")
  md <- obj[[]]
  md$cluster <- as.character(Idents(obj))

  cluster_broad <- aggregate(md[, broad_cols, drop = FALSE],
                             by = list(cluster = md$cluster), FUN = mean)
  cluster_broad$broad_lineage <- names(broad_modules)[
    max.col(as.matrix(cluster_broad[, broad_cols, drop = FALSE]), ties.method = "first")
  ]

  cluster_epi <- aggregate(md[, epi_cols, drop = FALSE],
                           by = list(cluster = md$cluster), FUN = mean)
  cluster_epi$epithelial_state <- names(epithelial_modules)[
    max.col(as.matrix(cluster_epi[, epi_cols, drop = FALSE]), ties.method = "first")
  ]

  obj$broad_lineage <- cluster_broad$broad_lineage[
    match(as.character(Idents(obj)), cluster_broad$cluster)
  ]
  obj$epithelial_state <- cluster_epi$epithelial_state[
    match(as.character(Idents(obj)), cluster_epi$cluster)
  ]
  obj$epithelial_state[obj$broad_lineage != "Epithelial"] <- "Non-epithelial"

  cluster_summary <- merge(cluster_broad, cluster_epi, by = "cluster")
  cluster_summary$day <- day_name
  write.csv(cluster_summary,
            file.path(out_dir, paste0(day_name, "_cluster_module_scores.csv")),
            row.names = FALSE)
  obj
}

day_objects <- Map(annotate_clusters, day_objects, names(day_objects))

marker_genes <- c(
  "TP63", "KRT5", "KRT14", "ITGA6", "IVL", "TGM1", "CLDN1",
  "SPRR1B", "ABCA12", "SPRR1A", "FLG", "KRT1", "KRT10", "DSG1"
)

summarize_markers <- function(obj, day_name) {
  genes <- intersect(marker_genes, rownames(obj))
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  meta <- obj[[]]
  groups <- interaction(meta$sample, meta$broad_lineage, meta$epithelial_state,
                        sep = "||", drop = TRUE)
  result <- lapply(levels(groups), function(gr) {
    cells <- which(groups == gr)
    parts <- strsplit(gr, "\\|\\|")[[1]]
    data.frame(
      day = day_name,
      sample = parts[1],
      broad_lineage = parts[2],
      epithelial_state = parts[3],
      gene = genes,
      n_cells = length(cells),
      average_log_normalized_expression = Matrix::rowMeans(expr[, cells, drop = FALSE]),
      percent_expressing = Matrix::rowMeans(expr[, cells, drop = FALSE] > 0) * 100,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

marker_summary <- do.call(rbind, Map(summarize_markers, day_objects, names(day_objects)))
write.csv(marker_summary,
          file.path(out_dir, "qPCR_candidate_scRNA_specificity_summary.csv"),
          row.names = FALSE)

cell_composition <- do.call(rbind, lapply(names(day_objects), function(day_name) {
  md <- day_objects[[day_name]][[]]
  tab <- as.data.frame(table(day = rep(day_name, nrow(md)), sample = md$sample,
                             broad_lineage = md$broad_lineage,
                             epithelial_state = md$epithelial_state),
                       stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, ]
  names(tab)[names(tab) == "Freq"] <- "n_cells"
  tab$sample_fraction <- ave(tab$n_cells, tab$sample, FUN = function(z) z / sum(z))
  tab
}))
write.csv(cell_composition, file.path(out_dir, "cell_composition.csv"), row.names = FALSE)

pdf(file.path(out_dir, "UMAP_annotations.pdf"), width = 11, height = 8.5)
for (day_name in names(day_objects)) {
  obj <- day_objects[[day_name]]
  print(DimPlot(obj, group.by = "broad_lineage", label = TRUE, repel = TRUE) +
          ggtitle(paste(day_name, "broad lineage")))
  print(DimPlot(obj, group.by = "epithelial_state", label = TRUE, repel = TRUE) +
          ggtitle(paste(day_name, "epithelial state")))
  print(DimPlot(obj, group.by = "sample") + ggtitle(paste(day_name, "sample")))
}
dev.off()

pdf(file.path(out_dir, "qPCR_candidate_dotplots.pdf"), width = 13, height = 7)
for (day_name in names(day_objects)) {
  obj <- day_objects[[day_name]]
  Idents(obj) <- factor(obj$epithelial_state,
                        levels = c("Basal", "Periderm", "Cycling", "Intermediate",
                                   "Terminal", "Non-epithelial"))
  print(DotPlot(obj, features = intersect(marker_genes, rownames(obj)), assay = "RNA") +
          RotatedAxis() + ggtitle(paste(day_name, "qPCR candidate expression")))
}
dev.off()

metadata_rows <- lapply(names(day_objects), function(day_name) {
  obj <- day_objects[[day_name]]
  umap <- Embeddings(obj, "umap")
  md <- obj[[]]
  data.frame(
    cell = rownames(md),
    day = day_name,
    sample = md$sample,
    cell_line = md$cell_line,
    cluster = as.character(Idents(obj)),
    broad_lineage = md$broad_lineage,
    epithelial_state = md$epithelial_state,
    nFeature_RNA = md$nFeature_RNA,
    nCount_RNA = md$nCount_RNA,
    percent_mt = md$percent.mt,
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    stringsAsFactors = FALSE
  )
})
metadata_all <- do.call(rbind, metadata_rows)
gz <- gzfile(file.path(out_dir, "cell_metadata_and_umap.csv.gz"), "wt")
write.csv(metadata_all, gz, row.names = FALSE)
close(gz)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Finished. Results: ", out_dir)
