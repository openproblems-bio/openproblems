glibrary(SingleCellExperiment)
library(Matrix)
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
library(rsvd)

# save the matrix of the obsm variable train

train_matrix <- reducedDim(sce, "train")

# this is the format used in the alraSeurat2.R example on KlugerLab/ALRA/ repo
reducedDim(sce, "train") <- Matrix(t(alra(t(as.matrix(train_matrix)))[[3]]))

# return

sce
