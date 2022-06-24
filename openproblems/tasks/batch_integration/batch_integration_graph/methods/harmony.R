library(SingleCellExperiment)
library(harmony)
library(Seurat)
seu <- as.Seurat(sce, data = NULL)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = rownames(seu@assays$originalexp), npcs = n_pca)
seu <- RunHarmony(seu,
  batch,
  max.iter.harmony = max_iter_harmony, max.iter.cluster = max_iter_cluster
)
sce <- as.SingleCellExperiment(seu)
sce
