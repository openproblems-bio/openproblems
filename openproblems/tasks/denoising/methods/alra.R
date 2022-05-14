# we explicitly cast to as.matrix, so unless as.matrix fails we don't need
# ALRA to check the class type of train_matrix. therefore, use fork of
# ALRA that doesn't check type:
source("https://raw.githubusercontent.com/wes-lewis/ALRA/master/alra.R")

# fetch the matrix of the obsm variable train_norm
train_matrix <- t(as.matrix(sce))

# sometimes ALRA randomly fails with a decomposition error; retry
alra_output <- NULL
attempt <- 0
while (is.null(alra_output)) {
  if (attempt < 4) {
    try(
      alra_output <- alra(train_matrix)
    )
  } else {
    # final attempt
    alra_output <- alra(train_matrix)
  }
  attempt <- attempt + 1
}

alra_mat <- alra_output$A_norm_rank_k_cor_sc

out <- t(as.matrix(alra_mat))

# return
out
