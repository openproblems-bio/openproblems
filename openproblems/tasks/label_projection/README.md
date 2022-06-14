# Label Projection

## The task

A major challenge for integrating single cell datasets is creating matching cell type annotations for each cell. One of the most common strategies for annotating cell types is referred to as ["cluster-then-annotate"](https://www.nature.com/articles/s41576-018-0088-9) whereby cells are aggregated into clusters based on feature similarity and then manually characterized based on differential gene expression or previously identified marker genes. Recently, methods have emerged to build on this strategy and annotate cells using [known marker genes](https://www.nature.com/articles/s41592-019-0535-3). However, these strategies pose a difficulty for integrating atlas-scale datasets as the particular annotations may not match.

To ensure that the cell type labels in newly generated datasets match existing reference datasets, some methods align cells to a previously annotated [reference dataset](https://academic.oup.com/bioinformatics/article/35/22/4688/54802990) and then _project_ labels from the reference to the new dataset.

Here, we compare methods for annotation based on a reference dataset. The datasets consist of two or more samples of single cell profiles that have been manually annotated with matching labels. These datasets are then split into training and test batches, and the task of each method is to train a cell type classifer on the training set and project those labels onto the test set.

## The metrics
Metrics for label projection aim to characterize how well each classifer correctly assigns cell type labels to cells in the test set.

* **Accuracy**: Average number of correctly applied labels.
* **F1 score**: The [F1 score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html) can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.
* **Micro F1 score**: TODO

## API

Datasets should contain the following attributes:

* `adata.obs["labels"]` with ground truth celltype labels,
* `adata.obs["batch"]` with information of batches in the data, and
* `adata.obs["is_train"]` with a train vs. test split

It should be noted that datasets may only contain a single batch, or not contain discriminative batch information.

Methods should assign output celltype labels to `adata.obs['labels_pred']` using only the labels from the training data.

Note that the true labels are contained in `adata['labels']`.

Metrics can compare `adata['labels']` to `adata.obs['labels_pred']` using only the labels from the test data.
