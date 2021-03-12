# Label Projection

Label projection refers to the automatic identification of cell identity labels in a test or query dataset based on a reference datasets that typically contains curated, manually labeled cells. For a review, see [Abdelaal et al. (2019)](https://doi.org/10.1186/s13059-019-1795-z).

## API

Datasets should contain the following attributes:

* `adata.obs["labels"]` with ground truth celltype labels,
* `adata.obs["batch"]` with information of batches in the data, and
* `adata.obs["is_train"]` with a train vs. test split

Methods should assign output celltype labels to `adata.obs['labels_pred']` using only the labels from the training data.

Note that the true labels are contained in `adata['labels']`.

Metrics can compare `adata['labels']` to `adata.obs['labels_pred']` using only the labels from the test data.
