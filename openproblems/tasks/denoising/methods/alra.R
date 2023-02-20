library(SingleCellExperiment)
# we explicitly cast to as.matrix, so unless as.matrix fails we don't need
# ALRA to check the class type of train_matrix. therefore, use fork of
# ALRA that doesn't check type:
source("https://raw.githubusercontent.com/wes-lewis/ALRA/master/alra.R")

# fetch the matrix of the obsm variable train_norm
train_matrix <- reducedDim(sce, "train_norm")
train_matrix <- t(as.matrix(train_matrix))

alra_output <- alra(train_matrix)

alra_mat <- alra_output$A_norm_rank_k_cor_sc

out <- t(as.matrix(alra_mat))

# return
out
