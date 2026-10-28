source("R/pipeline.R")
suppressPackageStartupMessages({ library(HDCytoData); library(readxl) })
options(timeout = 600)

dir.create("data/PBMC8", recursive = TRUE, showWarnings = FALSE)

## 1. metadata / panel 表（Zenodo）
zen <- "https://zenodo.org/records/10039274/files"
for (f in c("PBMC8_metadata.xlsx", "PBMC8_panel_v3.xlsx")) {
  dest <- file.path("data/PBMC8", f)
  if (!file.exists(dest))
    download.file(paste0(zen, "/", f, "?download=1"), dest, mode = "wb")
}
md    <- as.data.frame(read_excel("data/PBMC8/PBMC8_metadata.xlsx"))
panel <- as.data.frame(read_excel("data/PBMC8/PBMC8_panel_v3.xlsx"))
print(md); print(head(panel))

## 2. 表达数据（HDCytoData / ExperimentHub，首次需几分钟，之后有缓存）
fs <- Bodenmiller_BCR_XL_flowSet()

## 3. 对齐检查：不对就停下来，别带着错的映射往下跑
print(flowCore::sampleNames(fs))
stopifnot(all(panel$fcs_colname %in% colnames(fs)))
if (!all(md$file_name %in% flowCore::sampleNames(fs)))
  stop("md$file_name 和 flowSet 的样本名对不上，把上面两个 print 的输出贴给我")

## 4. 因子顺序：Ref 在前，对比方向才是 BCRXL vs Ref
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])

## 5. 构建 SCE 并跑流程
sce <- prepData(fs, panel, md, features = panel$fcs_colname)
print(table(sce$sample_id))
run_pipeline(sce, tag = "pbmc8_full")
