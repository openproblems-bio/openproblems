# Label Projection

Label projection refers to the automatic identification of cell types in a test set based on a curated set of manually labeled cells. For a review, see [Abdelaal et al. (2019)](https://doi.org/10.1186/s13059-019-1795-z).

## API

Datasets should contain the following attributes:

* `adata.obs["labels"]` (ground truth celltype labels)
* `adata.obs["is_train"]` (train vs. test boolean)

Methods should assign celltype labels to `adata.obs['labels_pred']` using only the labels from the training data. The true labels are contained in `adata['labels']`.

Metrics should compare `adata['labels']` to `adata.obs['labels_pred']` using on the labels from the test data.
