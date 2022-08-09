# Universal Data Loader

Download data from URL and ensure that matrices are stored as expected.

## API
Downloaded datasets contain:

* `adata.X`: raw counts
* `adata.layers['lognorm']`: log-normalized counts, only if available
