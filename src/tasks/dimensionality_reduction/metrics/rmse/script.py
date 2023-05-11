import anndata as ad
import numpy as np
import sklearn.decomposition
import scipy.optimize
import scipy.spatial
from sklearn.metrics import pairwise_distances
import umap
import umap.spectral

## VIASH START
par = {
    "input_reduced": "resources_test/dimensionality_reduction/pancreas/reduced.h5ad",
    "input_test": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "score.h5ad",
}
## VIASH END

def _rmse(X, X_emb):
    high_dimensional_distance_vector = scipy.spatial.distance.pdist(X)
    low_dimensional_distance_vector = scipy.spatial.distance.pdist(X_emb)
    _, rmse = scipy.optimize.nnls(
        low_dimensional_distance_vector[:, None], high_dimensional_distance_vector
    )
    return rmse

print("Load data", flush=True)
input_test = ad.read_h5ad(par["input_test"])
input_reduced = ad.read_h5ad(par["input_reduced"])

high_dim = input_test.layers["normalized"]
X_emb = input_reduced.obsm["X_emb"]

print("Compute NNLS residual after SVD", flush=True)
n_svd = 200
svd_emb = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(high_dim)
rmse = _rmse(svd_emb, X_emb)

print("Compute NLSS residual after spectral embedding", flush=True)
n_comps = min(200, min(input_test.shape) - 2)
umap_graph = umap.UMAP(transform_mode="graph").fit_transform(high_dim)
spectral_emb = umap.spectral.spectral_layout(
    high_dim, umap_graph, n_comps, random_state=np.random.default_rng()
)
rmse_spectral = _rmse(spectral_emb, X_emb)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_test.uns["dataset_id"],
        "normalization_id": input_test.uns["normalization_id"],
        "method_id": input_reduced.uns["method_id"],
        "metric_ids": [ "rmse", "rmse_spectral" ],
        "metric_values": [ rmse, rmse_spectral ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")