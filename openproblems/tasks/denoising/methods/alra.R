library(SingleCellExperiment)
library(Matrix)
source("https://raw.githubusercontent.com/KlugerLab/ALRA/v1.0.0/alra.R")
library(rsvd)

# save the matrix of the obsm variable train

train_matrix <- t(as.matrix(reducedDim(sce, "train")))

alra_output <- alra(as.matrix(train_matrix))
alra_mat <- alra_output$A_norm_rank_k_cor_sc

reducedDim(sce, "train") <- as(t(alra_mat), "dgCMatrix")

# return
sce
