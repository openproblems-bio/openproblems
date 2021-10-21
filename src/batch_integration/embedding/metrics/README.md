# Evaluation of Batch Integration with Embedding Output

Metrics on embedding output include:

* Average silhouette width on batches ASW_batch
* Average silhouette width on labels ASW_label
* Cell cycle conservation
* Principle component regression PCR

## API

All datasets should contain the following attributes:

* `adata.uns['name']`: name of the dataset
* `adata.obs['batch']`: the batch covariate
* `adata.obs['label']`: the cell identity label
* `adata.obsm['X_pca']`: the PCA embedding before integration
* `adata.obsm['X_emb']`: the embedding after integration

Metrics compare:

* `adata.obsm['X_emb']` to `adata.obsm['X_pca']`
* `adata.obsm['X_emb']` to `adata.obs['label']`
* `adata.obsm['X_emb']` to `adata.obs['batch']`
