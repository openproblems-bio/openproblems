+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = ""
subtitle = ""

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"

[design.background]
  # Text color (true=light or false=dark).
  text_color_light = false
  color = "white"

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS.
 css_style = ""

 # CSS class.
 css_class = ""


+++
## Our mission
Our goal is to facilitate the development of novel computational methods to address open problems in the single-cell field. We are focused on bridging the gap between experts in computer science and machine learning and the biological problems associated with the single cell data. We want to identify important problems, aggregate standardized datasets, and create a platform to benchmark novel methods against the current state of the art using a common set of test metrics.

## Who can get involved
We want to build a diverse and inclusive community to support the Open Problems. As such we welcome any individual who wants to get involved and agrees to follow our [Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). We are currently supported by the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/) and welcome participation from labs across the single cell and/or machine learning communities, and in particular labs already involved in the [Single Cell Biology Seed Networks](https://chanzuckerberg.com/science/programs-resources/single-cell-biology/seednetworks/).

## How the Open Problems are structured

We have broken down the development of Single Cell Open Problems into **tasks**. A task is a specific quantifiable problem that addresses an open problem in the single-cell field. An example of a task is [Multimodal Data Integration](results/#multimodal_data_integration) in which the goal is to align single-cell measurements of different -omics modalities that will enable us to build increasingly complex characterizations of cell types and states.
#We evaluate these methods by using datasets where multimodal data is measured on the exact same cells (e.g. joint single cell RNA and ATAC profiling) to assess which methods correctly matched these measurements without using cell barcodes.

Each task is composed of three components:

* **Datasets** - a group of high quality curated datasets matched to a task
* **Metrics** - a set of quantitative measures that are used to rank methods
* **Methods** - algorithms contributed by the community to perform the task

All of the code for these tasks are hosted in an open source [GitHub repository](https://github.com/openproblems-bio/openproblems).

## Adding a dataset

### Data downloaders vs data loaders

Datasets are collections of single cell measurements that can be used for benchmarking a task. To add a dataset, you need to set up two components:

1. A data downloader function in `openproblems/data`
2. A task-specific data loader function in `openproblems/task/<task_name>/datasets`

The role of the data downloader is to grab the data from a public repository, perform any necessary preprocessing, and return an [`AnnData`](https://github.com/theislab/anndata/) object.

The API of a **data downloader** is

```
function dataset(bool test=False) -> AnnData adata
```

If `test` is True, then the method should load the full dataset, but only return a small version of the same data (preferably <200 cells and <500 genes) that can be used for testing purposes. The loaded AnnData objects are then used to evaluate various methods.

Next, we need a task-specific **data loader** that loads the data in a way that's formatted correctly for a given task. The specific data format for each task can be found in the `README.md` file in each `openproblems/tasks/<task_name>` directory. Generally speaking, `adata.X` should contain UMI counts (or equivalent). For example, the label projection task has the following requirements:


> ## API
>
> Datasets should contain the following attributes:
>
> * `adata.obs["labels"]` (ground truth celltype labels)
> * `adata.obs["is_train"]` (train vs. test boolean)

Note, we may be able to prepare a single dataset multiple ways for a single task. For example, in the zebrafish dataset that comprises single cell profiles from two different labs, we can create a train/test split based on the lab or by randomly splitting cells regardless of where they were measured.

## Adding a metric

Metrics are used to compare the output of each method and are task-specific. You can find a thorough discussion of model evaluation metrics in the `sklearn` [User Guide](https://scikit-learn.org/stable/modules/model_evaluation.html). The API of a metric is as follows:

```
function metric(adata) -> float
```

We encourage developers to submit a variety of metrics for each task since each method for model evaluation has specific biases. For example, [a recent comparison](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full) of dataset integration methods used 14 different evaluation metrics to compare methods.

## Adding a method

### Introduction to methods

Methods are the backbone of Single Cell Open Problems, and we hope that this is where most of the development will occur. Like metrics, methods are task-specific. The exact API of each method can be found in each `openproblems/task/<task_name>/methods/` directory. For example, the label projection task has the following API

> Methods should assign celltype labels to `adata.obs['labels_pred']` using only the labels from the training data. The true labels are contained in `adata['labels']`.

```
function _labelprojection(adata) -> adata.obs['labels_pred']
```

The `_` precedes the method name because it is not intended to be called directly during benchmarking. Instead, this method will be combined with preprocessing functions to create a full method as described in the following section.

### Handling preprocessing

Preprocessing is a major factor affecting the output of a single cell analysis pipeline, yet there is [little consensus](https://doi.org/10.1186/s13059-020-02136-7) on the optimal set of preprocessing steps for any given single cell task. Our approach is to provide multiple preprocessing options that can be easily combined with any given method.

Functions for normalization (accounting for varying UMIs per cell) and transformation (scaling for differences in average detection of each gene) can be found in [`openproblems/tools/normalize.py`](https://github.com/openproblems-bio/openproblems/blob/main/openproblems/tools/normalize.py).

We currently provide three flavors of normalization:
* `log_cpm` - log-transformed, counts-per-million normalized
* `sqrt_cpm` - sqrt-transformed, counts-per-million normalized
* `log_scran_pooling` - log-transformed, [scran](https://doi.org/10.1186/s13059-016-0947-7) normalized

To define a full method, we need to combine the base method above with a preprocessing function. For example, in the [`logistic_regression.py`](https://github.com/openproblems-bio/openproblems/blob/main/openproblems/tasks/label_projection/methods/logistic_regression.py) script, we define a `_logistic_regression()` base function and then combine it with preprocessing steps in the `logistic_regression_log_cpm()` and `logistic_regression_scran()` functions.

### Pull requests trigger benchmarking

Our current infrastructure will evaluate the performance of a method once the code has been merged to the `main` branch. We encourage you to take advantage of the test versions of each dataset and ensure that a new method will run properly on each dataset. Next, the easiest way to evaluate the performance of the method against all the datasets assigned to a task is to submit a pull request on GitHub.

## Adding a task

To add a new task, you need to collect all three of the above components of a task:

* One or more **datasets**
* One or more **metrics**
* One or more **methods**

We'd love to see new tasks added to the framework, and our core group of developers can help get new tasks off the ground. We already have some proposed tasks in our [GitHub Issues tracker](https://github.com/openproblems-bio/openproblems/issues?q=is%3Aissue+is%3Aopen+label%3Atask). Join us on Slack to get started!


## Join us on Slack

If you have any questions or would like to get more involved, please join us on the [CZI Science Slack](https://join-cziscience-slack.herokuapp.com/). Once you've created an account, look for the `#open-problems-sca` channel.
