from ....tools.decorators import metric

import copy
import harmonypy
import numpy as np
import scanpy as sc
import sklearn


@metric(metric_name="Label transfer distance", maximize=False)
def label_transfer_dist(adata):
    # Create merge annData with imputed and real gene activity
    val = adata
    test = copy.deepcopy(adata)
    test.X = np.array(test.obsm["gene_score"].todense())
    adata_merge = val.concatenate(test, batch_categories=["val", "test"])
    # Run harmony
    sc.tl.pca(adata_merge)
    sc.external.pp.harmony_integrate(adata_merge, "batch")
    # Compute distance in batch-corrected PCA embedding
    dist = sklearn.metrics.pairwise.euclidean_distances(
        adata_merge[adata_merge.obs.batch == "val"].obsm["X_pca_harmony"],
        adata_merge[adata_merge.obs.batch == "test"].obsm["X_pca_harmony"],
    )
    # Mean distance between cell pairs
    result = np.mean(np.triu(dist))
    return result
