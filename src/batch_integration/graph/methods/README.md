# Run Graph Integration Methods

Viash component for running integration methods.

## API

### Input data formats

The components before integration must contain:

* `adata.uns['name']`: name of the dataset
* `adata.obs['batch']`: batch covariate
* `adata.X`: log-normalized expression

Whether a dataset is scaled or selected for highly-variable genes before integration is passed by parameters to the
script.

### Output data formats

* `adata.obsp['connectivities']`: Integrated graph connectivities
* `adata.obsp['distances']`: Integrated graph distances
* `adata.uns['hvg']`: Number of highly variable genes selected before integration (0 meaning that no HVG selection was performed)
* `adata.uns['scaled']`: Boolean entry whether scaling was performed before integration
