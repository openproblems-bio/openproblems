library(SingleCellExperiment)
library(liger)
library(Seurat)

# Only counts is converted to liger object. To pass our own normalized data,
# store it in the "counts" slot
sobj <- as.Seurat(sce, data = NULL)
sobj@assays$RNA <- sobj@assays$originalexp # nolint
sobj@assays$RNA@counts <- sobj@assays$RNA@data

# Create Liger object
lobj <- seuratToLiger(sobj,
  combined.seurat = T, meta.var = batch, renormalize = FALSE,
  remove.missing = FALSE
)

# We only pass nomarlized data, so store it as such
lobj@norm.data <- lobj@raw.data

# Assign hvgs
lobj@var.genes <- rownames(sobj@assays$RNA)

lobj <- scaleNotCenter(lobj, remove.missing = FALSE)

# Use tutorial coarse k suggests of 20.
lobj <- optimizeALS(lobj, k = k, thresh = thresh, nrep = nrep)

lobj <- quantileAlignSNF(lobj,
  resolution = 0.4,
  small.clust.thresh = 20
)

# Return embedding
lobj@H.norm
