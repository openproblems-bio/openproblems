cat("Load dependencies\n")
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_test_mod2 = "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  input_prediction = "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
  output = "openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad"
)
#/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/work/29/320fe1e10fcd323020345bcc8969c2/openproblems_bmmc_cite_mod2_dummy_mean_per_gene.correlation.output.h5ad
## VIASH END

cat("Reading solution file\n")
ad_sol <- anndata::read_h5ad(par$input_test_mod2)

cat("Reading prediction file\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)

cat("Check prediction format\n")
expect_equal(
  ad_sol$uns$dataset_id, ad_pred$uns$dataset_id,
  info = "Prediction and solution have differing dataset_ids"
)

expect_true(
  isTRUE(all.equal(dim(ad_sol), dim(ad_pred))),
  info = "Dataset and prediction anndata objects should have the same shape / dimensions."
)

cat("Computing correlation metrics\n")
# Wrangle data
tv <- ad_sol$layers[["normalized"]]
pv <- ad_pred$layers[["normalized"]]

# precompute sds
tv_sd2 <- proxyC::colSds(tv)
pv_sd2 <- proxyC::colSds(pv)
tv_sd1 <- proxyC::rowSds(tv)
pv_sd1 <- proxyC::rowSds(pv)

# Compute metrics
pearson_vec_1 <- diag(dynutils::calculate_similarity(tv, pv, method = "pearson", margin = 1, diag = TRUE, drop0 = TRUE))
spearman_vec_1 <- diag(dynutils::calculate_similarity(tv, pv, method = "spearman", margin = 1, diag = TRUE, drop0 = TRUE))

pearson_vec_1[tv_sd1 == 0 | pv_sd1 == 0] <- 0
spearman_vec_1[tv_sd1 == 0 | pv_sd1 == 0] <- 0
# pearson_vec_1[!is.finite(pearson_vec_1) | pearson_vec_1 > 10] <- 0
# spearman_vec_1[!is.finite(spearman_vec_1) | spearman_vec_1 > 10] <- 0

mean_pearson_per_cell <- mean(pearson_vec_1)
mean_spearman_per_cell <- mean(spearman_vec_1)

pearson_vec_2 <- diag(dynutils::calculate_similarity(tv, pv, method = "pearson", margin = 2, diag = TRUE, drop0 = TRUE))
spearman_vec_2 <- diag(dynutils::calculate_similarity(tv, pv, method = "spearman", margin = 2, diag = TRUE, drop0 = TRUE))

pearson_vec_2[tv_sd2 == 0 | pv_sd2 == 0] <- 0
spearman_vec_2[tv_sd2 == 0 | pv_sd2 == 0] <- 0
# pearson_vec_2[!is.finite(pearson_vec_2) | pearson_vec_2 > 10] <- 0
# spearman_vec_2[!is.finite(spearman_vec_2) | spearman_vec_2 > 10] <- 0

mean_pearson_per_gene <- mean(pearson_vec_2)
mean_spearman_per_gene <- mean(spearman_vec_2)

overall_pearson <- cor(as.vector(tv), as.vector(pv), method = "pearson")
overall_spearman <- cor(as.vector(tv), as.vector(pv), method = "spearman")

metric_ids <- c("mean_pearson_per_cell", "mean_spearman_per_cell", "mean_pearson_per_gene", "mean_spearman_per_gene", "overall_pearson", "overall_spearman")
metric_values <- c(mean_pearson_per_cell, mean_spearman_per_cell, mean_pearson_per_gene, mean_spearman_per_gene, overall_pearson, overall_spearman)

cat("Create output object\n")
out <- anndata::AnnData(
  obs = data.frame(row.names = rownames(ad_sol), pearson = pearson_vec_1, spearman = spearman_vec_1),
  var = data.frame(row.names = colnames(ad_sol), pearson = pearson_vec_2, spearman = spearman_vec_2),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = metric_ids,
    metric_values = metric_values
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
