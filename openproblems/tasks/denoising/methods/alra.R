glibrary(SingleCellExperiment)
library(Matrix)
source("https://raw.githubusercontent.com/wes-lewis/ALRA/master/alra.R")
library(rsvd)

# save the matrix of the obsm variable train

train_matrix <- t(as.matrix(reducedDim(sce, "train")))

# cast to as.matrix, so unless as.matrix fails we don't need ALRA to check
# the class type of train_matrix. therefore, use version of ALRA that
# doesn't check type: https://github.com/wes-lewis/ALRA/blob/master/alra.R
alra_output <- alra(as.matrix(train_matrix))
alra_mat <- alra_output$A_norm_rank_k_cor_sc

reducedDim(sce, "train") <- as(t(alra_mat), "dgCMatrix")

# return
sce
