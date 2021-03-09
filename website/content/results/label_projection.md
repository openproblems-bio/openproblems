+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 30  # Order that this section will appear.

title = "Label projection"
subtitle = "*Using cell labels from a reference dataset to annotate an unseen dataset*"

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"


  # Text color (true=light or false=dark).
  text_color_light = false

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS.
 css_style = ""

 # CSS class.
 css_class = ""

datafile = "/assets/results/test.csv"
+++
## The task

A major challenge for integrating single cell datasets is creating matching cell type annotations for each cell. One of the most common strategies for annotating cell types is referred to as ["cluster-then-annotate"](https://www.nature.com/articles/s41576-018-0088-9) whereby cells are aggregated into clusters based on feature similarity and then manually characterized based on differential gene expression or previously identified marker genes. Recently, methods have emerged to build on this strategy and annotate cells using [known marker genes](https://www.nature.com/articles/s41592-019-0535-3). However, these strategies pose a difficulty for integrating atlas-scale datasets as the particular annotations may not match.

To ensure that the cell type labels in newly generated datasets match existing reference datasets, some methods align cells to a previously annotated [reference dataset](https://academic.oup.com/bioinformatics/article/35/22/4688/54802990) and then _project_ labels from the reference to the new dataset.

Here, we compare methods for annotation based on a reference dataset. The datasets consist of two or more samples of single cell profiles that have been manually annotated with matching labels. These datasets are then split into training and test batches, and the task of each method is to train a cell type classifer on the training set and project those labels onto the test set.

## The metrics
Metrics for label projection aim to characterize how well each classifer correctly assigns cell type labels to cells in the test set.

* **Accuracy**: Average number of correctly applied labels.
* **F1 score**: The [F1 score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html) can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.

## The results

{{< task_tables task="label_projection" >}}
