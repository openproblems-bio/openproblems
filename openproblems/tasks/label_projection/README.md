# Label Projection

Here's a brief task description, maybe link to some seminal papers.

## API

Datasets should contain the following attributes:

* `adata.obs["labels"]` (ground truth celltype labels)
* `adata.obs["is_train"]` (train vs. test boolean)

Methods should assign celltype labels to `adata.obs['labels_pred']` using only the labels from the training data. The true labels are contained in `adata['labels']`.

Metrics should compare `adata['labels']` to `adata.obs['labels_pred']` using on the labels from the test data.
