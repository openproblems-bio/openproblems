#' Seuratv3 TransferData wrapper
#'
#' rpy2 + seurat causes R to hang indefinitely.
#' This hack forks R to a separate process.
#'
#' @param sce SingleCellExperiment
#' @param n_pcs int Number of PCA components
#' @param k_filter int How many neighbors (k) to use when filtering anchors.
#' Set to NA to turn off filtering.
#' @param k_score int How many neighbors (k) to use when scoring anchors
#' @param script_path character Path to seuratv3.R

saveRDS(
  list(sce = sce, n_pcs = n_pcs, k_filter = k_filter, k_score = k_score),
  "/tmp/openproblems_seurat_args.rds"
)
# clear memory
rm("sce")

exitcode <- system2("Rscript", script_path)
if (exitcode != 0) {
  stop(paste0("Rscript failed with error code ", as.character(exitcode)))
}

sce_sp <- readRDS("/tmp/openproblems_seurat_sce_out.rds")
sce_sp
