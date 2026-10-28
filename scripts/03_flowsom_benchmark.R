suppressPackageStartupMessages({
  library(HDCytoData); library(SummarizedExperiment); library(flowCore)
  library(FlowSOM); library(mclust); library(clue)
})
set.seed(1234)
dir.create("results/benchmark", recursive = TRUE, showWarnings = FALSE)

## 一对一匹配的 F1（匈牙利算法），未匹配到的参考群体记 0
hungarian_f1 <- function(ref, pred) {
  tab <- unclass(table(ref, pred))
  Fm  <- 2 * tab / outer(rowSums(tab), colSums(tab), "+")
  n   <- max(dim(Fm))
  P   <- matrix(0, n, n); P[seq_len(nrow(Fm)), seq_len(ncol(Fm))] <- Fm
  a   <- as.integer(clue::solve_LSAP(P, maximum = TRUE))
  mean(P[cbind(seq_len(nrow(Fm)), a[seq_len(nrow(Fm))])])
}

datasets <- list(Samusik_01 = Samusik_01_SE, Levine_32dim = Levine_32dim_SE)

res <- do.call(rbind, lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]()
  is_type <- colData(d)$marker_class == "type"
  x   <- as.matrix(asinh(assay(d, "exprs")[, is_type] / 5))
  ref <- as.character(rowData(d)$population_id)

  ## 有人工标签的细胞：排除 NA 和 unassigned
  labeled <- !is.na(ref) & !grepl("unassigned", ref, ignore.case = TRUE)
  k <- length(unique(ref[labeled]))

  ## 1. 对全部细胞聚类
  ff <- flowFrame(x); cols <- seq_len(ncol(x))
  f_k  <- FlowSOM(ff, colsToUse = cols, xdim = 10, ydim = 10, nClus = k,  seed = 1234)
  f_40 <- FlowSOM(ff, colsToUse = cols, xdim = 10, ydim = 10, nClus = 40, seed = 1234)
  pred_k   <- as.character(GetMetaclusters(f_k))
  pred_40  <- as.character(GetMetaclusters(f_40))
  pred_som <- as.character(GetClusters(f_k))   # 100 个 SOM 节点

  ## 2. 只在有标签的细胞上评估
  r <- ref[labeled]
  ev <- function(p) p[labeled]

  data.frame(
    dataset = nm,
    n_cells_total = nrow(x),
    n_cells_labeled = sum(labeled),
    pct_unassigned = round(100 * (1 - mean(labeled)), 1),
    n_ref_populations = k,
    setting = c("metaclusters = n_ref", "metaclusters = 40", "100 SOM nodes"),
    ARI = c(adjustedRandIndex(r, ev(pred_k)), NA, NA),
    mean_F1_hungarian = c(hungarian_f1(r, ev(pred_k)),
                          hungarian_f1(r, ev(pred_40)),
                          hungarian_f1(r, ev(pred_som)))
  )
}))
print(res)
write.csv(res, "results/benchmark/flowsom_vs_manual.csv", row.names = FALSE)
