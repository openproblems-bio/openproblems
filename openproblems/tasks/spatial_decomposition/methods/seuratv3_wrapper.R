#' Seuratv3 TransferData wrapper
#'
#' rpy2 + seurat causes R to hang indefinitely.
#' This hack forks R to a separate process.
#'
#' @param sce_sc SingleCellExperiment single-cell data
#' @param sce_sp SingleCellExperiment spatial data
#' @param n_pcs int Number of principal components
#' @param script_path character Path to seuratv3.R

saveRDS(
  list(sce_sc = sce_sc, sce_sp = sce_sp, n_pcs = n_pcs),
  "/tmp/openproblems_seurat_args.rds"
)
# clear memory
rm(c(sce_sc, sce_sp))

exitcode <- system2("Rscript", script_path)
if (exitcode != 0) {
  stop(paste0("Rscript failed with error code ", as.character(exitcode)))
}

sce_sp <- readRDS("/tmp/openproblems_seurat_sce_sp_out.rds")
sce_sp
