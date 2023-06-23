import anndata as ad
import numpy as np
import yaml


## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'functionality_name': 'foo',
    'config': 'bar',
}

## VIASH END

def _set_uns(adata):
    adata.uns["neighbors"] = adata.uns["knn"]
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    adata.uns["neighbors"]["distances_key"] = "distances"

def _randomize_features(X, partition=None):
    X_out = X.copy()
    if partition is None:
        partition = np.full(X.shape[0], 0)
    else:
        partition = np.asarray(partition)
    for partition_name in np.unique(partition):
        partition_idx = np.argwhere(partition == partition_name).flatten()
        X_out[partition_idx] = X[np.random.permutation(partition_idx)]
    return X_out

def _randomize_graph(adata, partition=None):
    distances, connectivities = (
        adata.obsp["knn_distances"],
        adata.obsp["knn_connectivities"],
    )
    new_idx = _randomize_features(np.arange(distances.shape[0]), partition=partition)
    adata.obsp["distances"] = distances[new_idx][:, new_idx]
    adata.obsp["connectivities"] = connectivities[new_idx][:, new_idx]
    _set_uns(adata)
    return adata

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["subtype"]

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])
input.X = input.layers["normalized"]

input.X = _randomize_features(input.X)
input.obsm["X_emb"] = _randomize_features(input.obsm["X_pca"])
input = _randomize_graph(input)
del input.X

print("Store outputs", flush=True)
input.uns['output_type'] = output_type
input.uns['method_id'] = meta['functionality_name']

input.write_h5ad(par['output'], compression='gzip')
