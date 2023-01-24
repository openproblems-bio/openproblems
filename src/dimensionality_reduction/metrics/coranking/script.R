library(anndata)
library(coRanking)

## VIASH START
par <- list(
  "input_reduced" = "resources_test/dimensionality_reduction/pancreas/reduced.h5ad",
  "input_test" = "resources_test/dimensionality_reduction/pancreas/test.h5ad",
  "output" = "score.h5ad"
)
## VIASH END

cat("Read anndata objects")
input_test <- anndata::read_h5ad(par[["input_test"]])
input_reduced <- anndata::read_h5ad(par[["input_reduced"]])

# get datasets
high_dim <- input_test$layers[["normalized"]]
X_emb <- input_reduced$obsm[["X_emb"]]

if (any(is.na(X_emb))) {
  continuity_at_k30 <-
    trustworthiness_at_k30 <-
    qnx_at_k30 <-
    lcmc_at_k30 <-
    qnx_auc <-
    qlocal <-
    qglobal <-
    0
} else {
  cat("Compute pairwise distances\n")
  # TODO: computing a square distance matrix is problematic for large datasets!
  # TODO: should we use a different distance metric for the high_dim?
  # TODO: or should we subset to the HVG?
  dist_highdim <- coRanking:::euclidean(as.matrix(high_dim))
  dist_emb <- coRanking:::euclidean(as.matrix(X_emb))

  cat("Compute ranking matrices\n")
  rmat_highdim <- rankmatrix(dist_highdim, input = "dist")
  rmat_emb <- rankmatrix(dist_emb, input = "dist")

  cat("Compute coranking matrix\n")
  corank <- coranking(rmat_highdim, rmat_emb, "rank")

  cat("Compute metrics\n")
  # Compute QNX. This is a curve indicating the percentage of points
  # that are mild in- and extrusions or keep their rank.
  qnx <- Q_NX(corank)

  # Calculate the local continuity meta-criterion from a co-ranking matrix.
  lcmc <- LCMC(corank)

  # the values of qnx are split into local and global values by kmax
  kmax <- which.max(lcmc)

  # check certain quality values at k=30
  k30 <- 30
  trustworthiness_at_k30 <- coRanking:::cm.M_T(corank, k30)
  continuity_at_k30 <- coRanking:::cm.M_C(corank, k30)
  qnx_at_k30 <- qnx[[k30]]
  lcmc_at_k30 <- lcmc[[k30]]

  # area under the QNX curve
  qnx_auc <- mean(qnx)

  # local quality measure
  qlocal <- mean(qnx[seq_len(kmax)])

  # global quality measure
  qglobal <- mean(qnx[-seq_len(kmax)])
}

cat("construct output AnnData\n")
output <- AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = input_test$uns[["dataset_id"]],
    normalization_id = input_test$uns[["normalization_id"]],
    method_id = input_reduced$uns[["method_id"]],
    metric_ids = c(
      "continuity_at_k30",
      "trustworthiness_at_k30",
      "qnx_at_k30",
      "lcmc_at_k30",
      "qnx_auc",
      "qlocal",
      "qglobal"
    ),
    metric_values = c(
      continuity_at_k30,
      trustworthiness_at_k30,
      qnx_at_k30,
      lcmc_at_k30,
      qnx_auc,
      qlocal,
      qglobal
    )
  )
)

cat("Write to file\n")
output$write_h5ad(par$output)
