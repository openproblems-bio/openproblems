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
Our goal is to facilitate the development of novel computational methods to address open problems in integrating the [Human Cell Atlas](https://www.humancellatlas.org/). We are focused on bridging the gap between experts in computer science and machine learning and the biological problems associated with the single cell data. We want to identify important problems, aggregate standardized datasets, and create a platform to benchmark novel methods against the current state of the art using a common set of test metrics.

## Who can get involved
We want to build a diverse and inclusive community the support the Open Problems. As such we welcome any individual who wants to get involved and agrees to follow our [Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). We are currently supported by the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/) and welcome participation from labs in the [Seed Networks for the Human Cell Atlas](https://chanzuckerberg.com/science/programs-resources/single-cell-biology/seednetworks/).

## How the Open Problems are structured

We have broken down the development of the Single Cell Open Problems into **tasks**. A task is a specific quantifiable problem that addresses an open problem in integrating the Human Cell Atlas. An example of a task is [Multimodal Data Integration](results/#multimodal_data_integration) in which the goal is to take a set of multimodal measurements of the same set of cells (e.g. joint single cell RNA and ATAC profiling) and identify which pairs of measurements were acquired from the same cell without using cell barcodes.

Each task is composed of three components:

* **Datasets** - a group of high quality curated datasets matched to a task
* **Metrics** - a set of quantitative measures that are used to rank methods
* **Methods** - algorithms contributed by the community to perform the task

All of the code for these tasks are hosted in an open source [GitHub repository](https://github.com/singlecellopenproblems/SingleCellOpenProblems). To understand how these components work together, let's examine the Multimodal Data Integration task.

## Multimodal Data Integration

### What's the task?
Several recently described technologies allow for simultaneous measurement of different aspects of cell state. For example, [sci-CAR](https://doi.org/10.1126/science.aau0730) jointly profiles RNA expression and chromatin accessibility on the same cell and [CITE-seq](https://doi.org/10.1038/nmeth.4380) measures surface protein abundance and RNA expression from each cell. However, these joint profiling methods have several tradeoffs compared to unimodal measurements. Joint methods can be more expensive or lower throughput or more noisy than measuring a single modality at a time. Therefore it is useful to develop methods that are capable of integrating measurements of the same biological system but obtained using different technologies. Here the goal is to match measurements acquired for the same cell without using the cell barcodes that link the two measurements.

### How's the GitHub organized?

To get started let's look at the structure of the [OpenProblems GitHub repository](https://github.com/singlecellopenproblems/SingleCellOpenProblems):

```python
openproblems/
├── data/ # functions to download raw datasets
│   ├── dummy.py
│   ├── __init__.py
│   ├── scicar/ # downloads sci-CAR data from the Gene Expression Omnibus (GEO)
│   │   ├── base.py
│   │   ├── cell_lines.py
│   │   ├── __init__.py
│   │   └── mouse_kidney.py
│   ├── utils.py
│   └── zebrafish.py
├── __init__.py
├── tasks/ # contains code specific to each task
│   ├── __init__.py
│   ├── label_projection/
│   └── `multimodal_data_integration`/ # all the code for the multimodal integration task
│       ├── checks.py
│       ├── datasets/ # functions to load task-specific versions of each dataset
│       ├── __init__.py
│       ├── methods/  # methods to perform multimodal integration
│       ├── metrics/  # metrics used to benchmark each task
│       └── README.md
├── test/ # unit tests for the module
├── tools/ # utility functions like normalization
├── utils.py
└── version.py
```

There are a few other things in the repository, such as the `website/` directory that contains the files the build the website you're reading right now, but for the purposes of getting involved, the `openproblems/` module is the most important.

The first thing that's important to note is that there is not necessarily a one-to-one correspondence between datasets and tasks. Intuitively, this makes sense if you consider that it may be useful to use a dataset both for data integration and visualization. Because of this, we separated the functions that download raw data from the code that prepares the data for a specific task.

Next, we look at the `tasks/` directory that contains the folder for the `multimodal_data_integration/` task.

```python
tasks/multimodal_data_integration/
├── checks.py
├── datasets/
│   ├── __init__.py
│   └── scicar.py
├── __init__.py
├── methods/
│   ├── cheat.py
│   ├── harmonic_alignment.py
│   ├── __init__.py
│   ├── mnn.py
│   └── procrustes.py
├── metrics/
│   ├── __init__.py
│   ├── knn_auc.py
│   └── mse.py
└── README.md
```

### What are the datasets?

To benchmark multimodal integration methods, we use joint-profiling datasets where we have a one-to-one correspondence between measurements of each data type. These matched datasets provide ground truth that can be used to determine how accurately a method can align measurements of the same biological system.
