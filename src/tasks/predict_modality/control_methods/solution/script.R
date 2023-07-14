cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_test_mod2 = "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  output = "output.h5ad"
)

meta <- list(
  functionality_name = "foo"
)
## VIASH END

cat("Reading h5ad files\n")
ad2_test <- anndata::read_h5ad(par$input_test_mod2)
ad2_test$uns[["method_id"]] <- meta$functionality_name

cat("Writing predictions to file\n")
zzz <- ad2_test$write_h5ad(par$output, compression = "gzip")
