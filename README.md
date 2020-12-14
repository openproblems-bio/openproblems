# SingleCellOpenProblems

![Travis CI Status](https://img.shields.io/travis/com/singlecellopenproblems/singlecellopenproblems/master?label=Travis%20CI)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/singlecellopenproblems/singlecellopenproblems/Run%20Tests/master?label=Github%20Actions)
[![Coverage Status](https://coveralls.io/repos/github/singlecellopenproblems/SingleCellOpenProblems/badge.svg?branch=master)](https://coveralls.io/github/singlecellopenproblems/SingleCellOpenProblems?branch=master)
[![Netlify Status](https://api.netlify.com/api/v1/badges/83b92388-53c7-4fef-9003-e14d94c6ac6f/deploy-status)](https://app.netlify.com/sites/openproblems/deploys)
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Style Guide: OpenStack](https://img.shields.io/badge/style%20guide-openstack-eb1a32.svg)](https://docs.openstack.org/hacking/latest/user/hacking.html#styleguide)

Formalizing and benchmarking open problems in single-cell genomics.

[**Visit the Open Problems Website**](https://openproblems.netlify.app/)

## Table of contents

- [The team](#the-team)
- [API](#api)
  * [Writing functions in R](#writing-functions-in-r)
  * [Adding package dependencies](#adding-package-dependencies)
- [Contributing](#contributing)
  * [Adding a new dataset](#adding-a-new-dataset)
  * [Adding a dataset / method / metric to a task](#adding-a-dataset---method---metric-to-a-task)
  * [Adding a new task](#adding-a-new-task)

<!-- Table of contents generated with [markdown-toc](http://ecotrust-canada.github.io/markdown-toc/) -->

## The team

**Core**:
* Scott Gigante (@scottgigante)  
* Daniel Burkhardt (@dburkhardt)
* Malte Luecken (@LuckyMD)  

**Task authors** (_alphabetically_):  
* Mohammad Lotfallahi (@M0hammadL) - Label projection task
* Qian Qin (@qinqian) - Predicting gene expression from chromatin accessibility  
* Michael Vinyard (@mvinyard) - Stress preservation in Dimensionality Reduction
* Florian Wagner (@flo-compbio) - Data denoising  

**Supervision** (_alphabetically_):  
* Smita Krishnaswamy, Yale  
* Fabian Theis, Helmholtz Munich  

**Chan Zuckerberg Initiative Support** (_alphabetically_):  
* Jonah Cool  
* Fiona Griffin

## API

Each task consists of datasets, methods, and metrics.

Datasets should take no arguments and return an AnnData object. If `test is True`, then the method should load the full dataset, but only return a small version of the same data (preferably <200 cells and <500 genes) for faster downstream analysis.

```
function dataset(bool test=False) -> AnnData adata
```

Methods should take an AnnData object and store the output in `adata.obs` according to the specification of the task.

```
function method(AnnData adata) -> AnnData adata
```

Metrics should take an AnnData object and return a `float`.

```
function metric(AnnData adata) -> float
```

Task-specific APIs are described in the README for each task.

* [Label Projection](openproblems/tasks/label_projection)
* [Multimodal Data Integration](openproblems/tasks/multimodal_data_integration)

### Writing functions in R

Metrics and methods can also be written in R, using [`scprep`'s `RFunction`](https://scprep.readthedocs.io/en/stable/reference.html#scprep.run.RFunction) class. AnnData Python objects are converted to and from `SingleCellExperiment` R objects using [`anndata2ri`](https://icb-anndata2ri.readthedocs-hosted.com/en/latest/). A simple method implemented in R could be written as follows:

```{python}
from ....tools.decorators import method
from ....tools.utils import check_version

import scprep

_pca = scprep.run.RFunction(
    setup="""
        library(SingleCellExperiment)
        library(stats)
    """,
    args="sce, n_pca=10",
    body="""
        reducedDim(sce, "pca") <- prcomp(t(assay(sce, "X")))[,1:n_pca]
        sce
    """,
)


@method(
    method_name="PCA",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_url="https://www.tandfonline.com/doi/abs/10.1080/14786440109462720",
    paper_year=1901,
    code_url="https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp",
    code_version=check_version("rpy2"),
    image="openproblems-r-base",
)
def pca(adata, n_pca=10):
    return _pca(adata, n_pca=n_pca)
```

See the [`anndata2ri` docs](https://icb-anndata2ri.readthedocs-hosted.com/en/latest/) for API details. For a more detailed example of how to use this functionality, see our implementation of [fastMNN batch correction](openproblems/tasks/multimodal_data_integration/methods/mnn.py).

### Adding package dependencies

If you are unable to write your method using our base dependencies, you may add to our existing Docker images, or create your own. The image you wish to use (if you are not using the base image) should be specified in the `image` keyword argument of the method/metric decorator. See the [Docker images README](docker/README.md) for details.

## Contributing

### Adding a new dataset

Datasets are loaded under `openproblems/data`. Each data loading function should download the appropriate dataset from a stable location (e.g. from Figshare) be decorated with `openproblems.data.utils.loader` in order to cache the result.

### Adding a dataset / method / metric to a task

To add a dataset, method, or metric to a task, simply create a new `.py` file corresponding to your proposed new functionality and import the main function in the corresponding `__init__.py`. E.g., to add a "F2" metric to the label projection task, we would create `openproblems/tasks/label_projection/metrics/f2.py` and add a line
```
from .f2 import f2
```
to [`openproblems/tasks/label_projection/metrics/__init__.py`](openproblems/tasks/label_projection/metrics/__init__.py).

For datasets in particular, these should be loaded using a `loader` function from `openproblems.data`, with only task-specific annotations added in the task-specific data file.

For methods and metrics, they should be decorated with the appropriate function in `openproblems.tools.decorators` to include metadata required for the evaluation and presentation of results.

Note that data is not normalized in the data loader; normalization should be performed as part of each method. For ease of use, we provide a collection of common normalization functions in [`openproblems.tools.normalize`](openproblems/tools/normalize.py).

### Adding a new task

The task directory structure is as follows

```
opensproblems/
  - tasks/
    - task_name/
      - README.md
      - __init__.py
      - checks.py
      - datasets/
        - __init__.py
        - dataset1.py
        - ...
      - methods/
        - __init__.py
        - method1.py
        - ...
      - metrics/
        - __init__.py
        - metric1.py
        - ...
```

`task_name/__init__.py` can be copied from an existing task.

`checks.py` should implement the following functions:

```
check_dataset(AnnData adata) -> bool # checks that a dataset fits the task-specific schema
check_method(AnnData adata) -> bool # checks that the output from a method fits the task-specific schema
```

For adding datasets, methods and metrics, see above.
