library(RcppDist)
library(anndata)

dest <- getwd()

# VIASH START
par <- list(
    "input_data" = "resources_test/spatially_variable_genes/mouse_brain_coronal_section1/dataset.h5ad",
    "output" = "output.h5ad",
    "n_iter" = 10
)
meta <- list(
    "functionality_name" = "BOOST-GP"
)
# VIASH END

cat("Load data\n")
adata <- anndata::read_h5ad(par$input_data)

setwd("/opt/BOOST-GP")
source("./R/boost.gp.R")

counts <- as.matrix(adata$layers[["counts"]])
colnames(counts) <- adata$var_names
rownames(counts) <- adata$obs_names
mode(counts) <- "integer"

loc <- as.data.frame(adata$obsm[["spatial"]])
rownames(loc) <- adata$obs_names
colnames(loc) <- c("x", "y")

cat("Run BOOST-GP\n")
df <- as.data.frame(boost.gp(Y = counts, loc = loc, iter = par$n_iter, burn = 5))

df$feature_id <- rownames(df)
df <- subset(df, select = c("feature_id", "PPI"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")

# save output
cat("Write output AnnData to file\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns[["dataset_id"]],
        "method_id" = meta[["functionality_name"]]
    )
)

zzz <- output$write_h5ad(paste0(dest, "/", par$output), compression = "gzip")
