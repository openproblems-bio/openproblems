# Packages
library(SingleCellExperiment)
library(scuttle)

# Input: sce
counts <- assay(sce, "X")

# Normalize
logNormCounts(sce, assay.type="X")

# TODO: update
colData(sce)$template_output <- 0

# Return
sce
