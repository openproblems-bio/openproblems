requireNamespace("scDesign3", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)
library(rlang)

# set random seed
set.seed(2024)

## VIASH START
par <- list(
  input = "resources_test/common/mouse_brain_coronal/dataset.h5ad",
  output = "dataset_sim.h5ad",
  gp_k = 50L,
  select_top_variable_genes = 50L
)
meta <- list(
  cpus = 30L
)
## VIASH END

cat("Read AnnData\n")
adata <- anndata::read_h5ad(par$input)

cat("Transform into SCE\n")
df_loc <- as.data.frame(adata$obsm[['spatial']])
colnames(df_loc) <- c("spatial1", "spatial2")
rownames(df_loc) <- adata$obs_names

ref_sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = Matrix::t(adata$layers[["counts"]])),
  colData = df_loc
)

ref_sce

# check the number of genes in reference object
n_genes <- dim(ref_sce)[1]

mu_formula <- paste0(
  "s(spatial1, spatial2, bs = 'gp', k = ", par$gp_k, ")"
)

if (n_genes > par$select_top_variable_genes) {
  cat("Select ", par$select_top_variable_genes, " genes among ", n_genes, " reference genes ", "\n", sep = "")

  cat("Transform into scDesign3 data format\n")
  ref_data <- scDesign3::construct_data(
    sce = ref_sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = NULL,
    spatial = c("spatial1", "spatial2"),
    other_covariates = NULL,
    corr_by = "1"
  )

  cat("Fit regression models for each feature\n")
  ref_marginal <- scDesign3::fit_marginal(
    data = ref_data,
    predictor = "gene",
    mu_formula = mu_formula,
    sigma_formula = "1",
    family_use = "nb",
    parallelization = "pbmcmapply",
    n_cores = 2L,
    usebam = FALSE,
    trace = TRUE
  )

  cat("Subset to the top variable genes\n")
  dev_explain <- sapply(ref_marginal, function(x) {
    if (length(x$fit) == 1 && is.na(x$fit)) {
      return(NA_real_)
    }
    summary(x$fit)$dev.expl
  })
  top_sel <- names(sort(dev_explain, decreasing = TRUE))[seq_len(par$select_top_variable_genes)]
} else {
  top_sel <- adata$var_names
}

ref_sce <- ref_sce[top_sel, ]
var_subset <- adata$var[top_sel, , drop = FALSE]

cat("Transform subset matrix into scDesign3 data format\n")
ref_data <- scDesign3::construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

cat("Fit expression of each gene with GP model\n")
ref_marginal <- scDesign3::fit_marginal(
  data = ref_data,
  predictor = "gene",
  mu_formula = mu_formula,
  sigma_formula = "1",
  family_use = "nb",
  parallelization = "pbmcmapply",
  n_cores = 2L,
  usebam = FALSE,
  trace = TRUE
)

cat("Fit a copula, obtain AIC and BIC\n")
ref_copula <- scDesign3::fit_copula(
  sce = ref_sce,
  assay_use = "counts",
  marginal_list = ref_marginal,
  family_use = "nb",
  copula = "gaussian",
  parallelization = "pbmcmapply",
  n_cores = 2L,
  input_data = ref_data$dat
)

cat("Extract out the estimated parameters\n")
ref_para <- scDesign3::extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat,
  parallelization = "pbmcmapply",
  n_cores = 2L
)

cat("Simulate the new count matrix\n")
# generate non-spatially variable mean values with shuffling
shuffle_idx <- sample(nrow(ref_para$mean_mat))
non_de_mat <- ref_para$mean_mat[shuffle_idx, ]

# simulate data with varied spatial variability
outputs <- lapply(seq(0, 1.0, 0.05), function(alpha){
  cat("Simulate data with alpha = ", alpha, "\n", sep = "")
  counts <- scDesign3::simu_new(
    sce = ref_sce,
    mean_mat = alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat,
    sigma_mat = ref_para$sigma_mat,
    zero_mat = ref_para$zero_mat,
    quantile_mat = NULL,
    copula_list = ref_copula$copula_list,
    n_cores = 5L,
    family_use = "nb",
    input_data = ref_data$dat,
    new_covariate = ref_data$newCovariate,
    important_feature = rep(TRUE, nrow(ref_sce)),
    filtered_gene = NULL
  )

  if ("feature_id" %in% names(var_subset)) {
    new_var <- data.frame(
      feature_id = paste0(var_subset$feature_id, "_", alpha),
      feature_name = paste0(var_subset$feature_name, "_", alpha),
      orig_feature_id = var_subset$feature_id,
      orig_feature_name = var_subset$feature_name,
      true_spatial_var_score = alpha
    )
    rownames(counts) <- new_var$feature_id
    rownames(new_var) <- new_var$feature_id
  } else {
    new_var <- data.frame(
      feature_id = paste0(var_subset$feature_name, "_", alpha),
      feature_name = paste0(var_subset$feature_name, "_", alpha),
      orig_feature_name = var_subset$feature_name,
      true_spatial_var_score = alpha
    )
    rownames(counts) <- new_var$feature_name
    rownames(new_var) <- new_var$feature_name
  }

  list(
    counts = Matrix::t(counts),
    var = new_var
  )
})

cat("Collecting final output\n", sep = "")
final_counts <- do.call(cbind, lapply(outputs, function(x) x$counts))
final_var <- do.call(rbind, lapply(outputs, function(x) x$var))
final_uns <- adata$uns[c("dataset_id", "dataset_name", "dataset_description", "dataset_summary", "dataset_url", "dataset_organism", "dataset_reference")]

output <- anndata::AnnData(
  layers = list(counts = final_counts),
  obs = adata$obs,
  var = final_var,
  obsm = adata$obsm,
  uns = final_uns
)

zzz <- output$write_h5ad(par$output, compression = "gzip")
