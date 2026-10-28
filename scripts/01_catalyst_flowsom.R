source("R/pipeline.R")
data(PBMC_fs, PBMC_panel, PBMC_md)
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
run_pipeline(sce, tag = "pbmc_demo")
