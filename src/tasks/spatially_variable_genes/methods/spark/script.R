suppressMessages(library(SPARK))
suppressMessages(library(anndata))

# VIASH START
par <- list(
    "input_data" = "resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad",
    "output" = "output.h5ad"
)
meta <- list(
    "functionality_name" = "SPARK",
    "cpus" = 4
)

# VIASH END

# load data
cat("Load data\n")
adata <- anndata::read_h5ad(par$input_data)
counts <- t(as.matrix(adata$layers[["counts"]]))
colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names
info <- as.data.frame(adata$obsm[["spatial"]])
rownames(info) <- colnames(counts)
colnames(info) <- c("x", "y")

# run SPARK
cat("Run SPARK\n")
if (!is.null(meta$cpus)) {
    n_cpus <- meta$cpus
} else {
    n_cpus <- 1
}

spark <- CreateSPARKObject(
    counts = counts, percentage = 0,
    min_total_counts = 0, location = info[, 1:2]
)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark,
    covariates = NULL,
    lib_size = spark@lib_size,
    num_core = n_cpus,
    verbose = FALSE
)

## Calculating pval
spark <- spark.test(spark,
    check_positive = T,
    verbose = F
)

df <- as.data.frame(spark@res_mtest)

df$feature_id <- rownames(df)

df <- subset(df, select = c("feature_id", "adjusted_pvalue"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")

# because SPARK only generates p-values, we here transform the values
# via -log10 to make sure a bigger score represents a higher spatial variation
df$pred_spatial_var_score <- -log10(df$pred_spatial_var_score)

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

anndata::write_h5ad(anndata = output, filename = par$output)
