## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input": "out_bash/multimodal_data_integration/methods/citeseq_cbmc_mnn.h5ad",
    "output": "out_bash/multimodal_data_integration/metrics/citeseq_cbmc_mnn_knn_auc.h5ad",
    "proportion_neighbors": 0.1,
    "n_svd": 100
}
## VIASH END

print("Importing libraries")
import anndata
import numpy as np
import sklearn.decomposition
import sklearn.neighbors

print("Reading adata file")
adata = anndata.read_h5ad(par["input"])

print("Checking parameters")
n_svd = min([par["n_svd"], min(adata.X.shape) - 1])
n_neighbors = int(np.ceil(par["proportion_neighbors"] * adata.X.shape[0]))

print("Performing PCA")
X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)

print("Compute KNN on PCA")
_, indices_true = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(X_pca)
    .kneighbors(X_pca)
)

print("Compute KNN on aligned matrix")
_, indices_pred = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(adata.obsm["aligned"])
    .kneighbors(adata.obsm["mode2_aligned"])
)

print("Check which neighbours match")
neighbors_match = np.zeros(n_neighbors, dtype=int)
for i in range(adata.shape[0]):
    _, pred_matches, true_matches = np.intersect1d(
        indices_pred[i], indices_true[i], return_indices=True
    )
    neighbors_match_idx = np.maximum(pred_matches, true_matches)
    neighbors_match += np.sum(
        np.arange(n_neighbors) >= neighbors_match_idx[:, None],
        axis=0,
    )

print("Compute area under neighbours match curve")
neighbors_match_curve = neighbors_match / (
    np.arange(1, n_neighbors + 1) * adata.shape[0]
)
area_under_curve = np.mean(neighbors_match_curve)

print("Store metic value")
adata.uns["metric_id"] = "knn_auc"
adata.uns["metric_value"] = area_under_curve

print("Writing adata to file")
adata.write_h5ad(par["output"], compression = "gzip")
