suppressPackageStartupMessages({
  library(CATALYST); library(diffcyt); library(flowCore)
  library(SingleCellExperiment); library(scater); library(ggplot2)
})

run_pipeline <- function(sce, tag, xdim = 10, ydim = 10, maxK = 20, seed = 1234) {
  set.seed(seed)
  figdir <- file.path("results", tag, "figures")
  tabdir <- file.path("results", tag, "tables")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

  fig <- function(p, name, w = 8, h = 6) {
    path <- file.path(figdir, name)
    if (inherits(p, "Heatmap")) {
      png(path, width = w, height = h, units = "in", res = 300)
      ComplexHeatmap::draw(p); invisible(dev.off())
    } else ggsave(path, p, width = w, height = h, dpi = 300)
  }

  ## QC
  fig(plotCounts(sce, group_by = "sample_id", color_by = "condition"), "01_counts.png")
  fig(pbMDS(sce, color_by = "condition"), "02_pbMDS.png")

  ## Clustering
  sce <- cluster(sce, features = "type", xdim = xdim, ydim = ydim,
                 maxK = maxK, seed = seed)
  k <- paste0("meta", maxK)
  fig(delta_area(sce), "03_delta_area.png")
  fig(plotExprHeatmap(sce, features = "type", by = "cluster_id", k = k,
                      bin_anno = TRUE), "04_cluster_heatmap.png", 10, 7)

  ## UMAP
  sce <- runDR(sce, "UMAP", cells = 500, features = "type")
  fig(plotDR(sce, "UMAP", color_by = k), "05_umap_clusters.png")
  fig(plotDR(sce, "UMAP", color_by = "condition"), "06_umap_condition.png")

  ## diffcyt (paired design: condition + patient)
  ei <- metadata(sce)$experiment_info
  design <- createDesignMatrix(ei, cols_design = c("condition", "patient_id"))
  print(design)
  contrast <- createContrast(c(0, 1, rep(0, ncol(design) - 2)))  # 第2列 = condition

  da <- diffcyt(sce, design = design, contrast = contrast,
                analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                clustering_to_use = k, verbose = FALSE)
  ds <- diffcyt(sce, design = design, contrast = contrast,
                analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                clustering_to_use = k, verbose = FALSE)

  write.csv(as.data.frame(rowData(da$res)), file.path(tabdir, "DA_results.csv"), row.names = FALSE)
  write.csv(as.data.frame(rowData(ds$res)), file.path(tabdir, "DS_results.csv"), row.names = FALSE)
  fig(plotDiffHeatmap(sce, rowData(da$res), all = TRUE, fdr = 0.05), "07_DA_heatmap.png")
  fig(plotDiffHeatmap(sce, rowData(ds$res), all = TRUE, fdr = 0.05), "08_DS_heatmap.png", 10, 8)

  invisible(sce)
}
