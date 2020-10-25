## Dimensional Reduction

Dimensional reduction is a key challenge in single-cell data representation.

## API

Datasets should contain the following attributes:

* `adata.obs["labels"]` (cell/observation labels)
* `adata.var["labels"]` (gene name feature labels)

**Methods** should assign dimensionally-reduced embedding coordinates to `adata.obsm['dimensional_reduction_method']`

**Metrics** should calculate the quality or "goodness of fit" of a dimensional reduction **method**. These matrices should be assigned to`adata.uns['metric']`.

