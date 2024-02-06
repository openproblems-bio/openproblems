cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod1.h5ad",
  input_mod2 = "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod2.h5ad",
  output_train_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod1.h5ad",
  output_train_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod2.h5ad",
  output_test_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/test_mod1.h5ad",
  output_test_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/test_mod2.h5ad",
  swap = TRUE,
  seed = 1L
)
# par <- list(
#   input_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/output_mod1.h5ad",
#   input_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/output_atac.h5ad",
#   output_train_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/train_mod1.h5ad",
#   output_train_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/train_mod2.h5ad",
#   output_test_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/test_mod1.h5ad",
#   output_test_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/test_mod2.h5ad",
#   swap = TRUE,
#   seed = 1L
# )
## VIASH END

cat("Using seed ", par$seed, "\n", sep = "")
set.seed(par$seed)

cat("Reading input data\n")
ad1 <- anndata::read_h5ad(if (!par$swap) par$input_mod1 else par$input_mod2)
ad2 <- anndata::read_h5ad(if (!par$swap) par$input_mod2 else par$input_mod1)

# figure out modality types
ad1_mod <- unique(ad1$var[["feature_types"]])
ad2_mod <- unique(ad2$var[["feature_types"]])

# determine new dataset id
new_dataset_id <- paste0(ad1$uns[["dataset_id"]], "_", tolower(ad1_mod), "2", tolower(ad2_mod))

# determine new uns
uns_vars <- c("dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism", "normalization_id")
ad1_uns <- ad2_uns <- ad1$uns[uns_vars]
ad1_uns$modality <- ad1_mod
ad2_uns$modality <- ad2_mod

# determine new obsm
ad1_obsm <- ad2_obsm <- list()

# determine new varm
ad1_var <- ad1$var[, intersect(colnames(ad1$var), c("gene_ids")), drop = FALSE]
ad2_var <- ad2$var[, intersect(colnames(ad2$var), c("gene_ids")), drop = FALSE]

if (ad1_mod == "ATAC") {
  # binarize features
  ad1$layers[["normalized"]]@x <- (ad1$layers[["normalized"]]@x > 0) + 0

  # copy gene activity in new object
  ad1_uns$gene_activity_var_names <- ad1$uns$gene_activity_var_names
  ad1_obsm$gene_activity <- as(ad1$obsm$gene_activity, "CsparseMatrix")
}

if (ad2_mod == "ATAC") {
  # subset to make the task computationally feasible
  if (ncol(ad2) > 10000) {
    poss_ix <- which(Matrix::colSums(ad2$layers[["normalized"]]) > 0)
    sel_ix <- sort(sample(poss_ix, 10000))
    ad2 <- ad2[, sel_ix]$copy()
    ad2_var <- ad2_var[sel_ix, , drop = FALSE]
  }

  # binarize features
  ad2$layers[["normalized"]]@x <- (ad2$layers[["normalized"]]@x > 0) + 0

  # copy gene activity in new object
  ad2_uns$gene_activity_var_names <- ad2$uns$gene_activity_var_names
  ad2_obsm$gene_activity <- as(ad2$obsm$gene_activity, "CsparseMatrix")
}

cat("Creating train/test split\n")
is_train <- which(ad1$obs[["is_train"]] == "train")
is_test <- which(!ad1$obs[["is_train"]] == "train")

# sample cells
if (length(is_test) > 1000) {
  ct <- as.character(ad1$obs[["cell_type"]][is_test])
  ct_tab <- table(ct)
  ct_freq <- setNames(as.vector(ct_tab) / sum(ct_tab), names(ct_tab))
  is_test <- sample(is_test, 1000, prob = sqrt(1 / ct_freq[ct]))
}

train_obs <- ad1$obs[is_train, intersect(colnames(ad1$obs), c("batch", "size_factors")), drop = FALSE]
test_obs <- ad1$obs[is_test, intersect(colnames(ad1$obs), c("batch", "size_factors")), drop = FALSE]
subset_mats <- function(li, obs_filt) {
  out <- list()
  for (n in names(li)) {
    out[[n]] <- li[[n]][obs_filt, , drop = FALSE]
  }
  out
}

cat("Create train objects\n")
output_train_mod1 <- anndata::AnnData(
  layers = subset_mats(list(counts = ad1$layers[["counts"]], normalized = ad1$layers[["normalized"]]), is_train),
  obsm = subset_mats(ad1_obsm, is_train),
  obs = train_obs,
  var = ad1_var,
  uns = ad1_uns
)
output_train_mod2 <- anndata::AnnData(
  layers = subset_mats(list(counts = ad2$layers[["counts"]], normalized = ad2$layers[["normalized"]]), is_train),
  obsm = subset_mats(ad2_obsm, is_train),
  obs = train_obs,
  var = ad2_var,
  uns = ad2_uns
)

cat("Create test objects\n")
output_test_mod1 <- anndata::AnnData(
  layers = subset_mats(list(counts = ad1$layers[["counts"]], normalized = ad1$layers[["normalized"]]), is_test),
  obsm = subset_mats(ad1_obsm, is_test),
  obs = test_obs,
  var = ad1_var,
  uns = ad1_uns
)
output_test_mod2 <- anndata::AnnData(
  layers = subset_mats(list(counts = ad2$layers[["counts"]], normalized = ad2$layers[["normalized"]]), is_test),
  obsm = subset_mats(ad2_obsm, is_test),
  obs = test_obs,
  var = ad2_var,
  uns = ad2_uns
)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")
