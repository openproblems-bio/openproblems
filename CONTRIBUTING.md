
# Contributing to Open Problems for Single Cell Analysis

There are many ways to contribute to `singlecellopenproblems`, with the most common ones
being contribution of datasets, methods, metrics, or new tasks. Improving the
documentation is no less important than improving the library itself. If you
find a typo in the documentation, or have made improvements, do not hesitate to
submit a GitHub pull request.

But there are many other ways to help. In particular answering queries on the
[issue tracker](https://github.com/KrishnaswamyLab/singlecellopenproblems/issues),
investigating bugs, and [reviewing other developers' pull requests](https://github.com/KrishnaswamyLab/singlecellopenproblems/pulls)
are very valuable contributions that decrease the burden on the project
maintainers.

Another way to contribute is to report issues you're facing, and give a "thumbs
up" on issues that others reported and that are relevant to you. It also helps
us if you spread the word: reference the project from your blog and articles,
link to it from your website, or simply star it in GitHub to say "I use it".

### Table of Contents

* [Submitting New Features](#submitting-new-features)
* [API](#api)
  + [Writing functions in R](#writing-functions-in-r)
  + [Adding package dependencies](#adding-package-dependencies)
  + [Adding a new dataset](#adding-a-new-dataset)
  + [Adding a dataset / method / metric to a task](#adding-a-dataset---method---metric-to-a-task)
  + [Adding a new task](#adding-a-new-task)
* [Code Style and Testing](#code-style-and-testing)
* [Code of Conduct](#code-of-conduct)
* [Attribution](#attribution)

<!-- Table of contents generated with [markdown-toc](http://ecotrust-canada.github.io/markdown-toc/) -->

## Submitting New Features

To submit new features to Open Problems for Single Cell Analysis, follow the steps below:

* Fork https://github.com/singlecellopenproblems/SingleCellOpenProblems
* Create repository secrets at [https://github.com/<username>/SingleCellOpenProblems/settings/secrets](https://github.com/<username>/SingleCellOpenProblems/settings/secrets):

*AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY are included in your AWS login details. If you do not have these, please contact us at [singlecellopenproblems@protonmail.com](mailto:singlecellopenproblems@protonmail.com).*

*TOWER_ACCESS_KEY (optional): log in with GitHub to https://tower.nf and create a token at https://tower.nf/tokens.*

* Enable workflows at [https://github.com/<username>/SingleCellOpenProblems/actions](https://github.com/<username>/SingleCellOpenProblems/actions)
* Set up your git repository:

```
git clone git@github.com:<username>/SingleCellOpenProblems.git
cd SingleCellOpenProblems
git remote add base git@github.com:singlecellopenproblems/SingleCellOpenProblems.git
git branch --set-upstream-to base/master
git pull
# IMPORTANT: choose a new branch name, e.g.
git checkout -b task/new_task_name # or metric/new_metric_name, etc
git push -u origin task/new_task_name
```

* Allow tests to pass on your new branch before pushing changes (as this will allow GitHub Actions to cache the workflow setup, which speeds up testing.)

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

Metrics and methods can also be written in R. AnnData Python objects are converted to and from `SingleCellExperiment` R objects using [`anndata2ri`](https://icb-anndata2ri.readthedocs-hosted.com/en/latest/). R methods should be written in a `.R` file which assumes the existence of a `SingleCellExperiment` object called `sce`, and should return the same object. A simple method implemented in R could be written as follows:

```{R}
### tasks/<task_name>/methods/pca.R
# Dependencies
library(SingleCellExperiment)
library(stats)

# Method body
n_pca <- 10
counts <- t(assay(sce, "X"))
reducedDim(sce, "pca") <- prcomp(counts)[,1:n_pca]

# Return
sce
```

```{python}
### tasks/<task_name>/methods/pca.py
from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

_pca = r_function("pca.R")


@method(
    method_name="PCA",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_url="https://www.tandfonline.com/doi/abs/10.1080/14786440109462720",
    paper_year=1901,
    code_url="https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp",
    code_version=check_version("rpy2"),
    image="openproblems-r-base",
)
def pca(adata):
    return _pca(adata)
```

See the [`anndata2ri` docs](https://icb-anndata2ri.readthedocs-hosted.com/en/latest/) for API details. For a more detailed example of how to use this functionality, see our implementation of fastMNN batch correction ([mnn.R](openproblems/tasks/multimodal_data_integration/methods/mnn.R), [mnn.py](openproblems/tasks/multimodal_data_integration/methods/mnn.py)).

### Adding package dependencies

If you are unable to write your method using our base dependencies, you may add to our existing Docker images, or create your own. The image you wish to use (if you are not using the base image) should be specified in the `image` keyword argument of the method/metric decorator. See the [Docker images README](docker/README.md) for details.

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
      - api.py
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

`api.py` should implement the following functions:

```
check_dataset(AnnData adata) -> bool # checks that a dataset fits the task-specific schema
check_method(AnnData adata) -> bool # checks that the output from a method fits the task-specific schema
sample_dataset() -> AnnData adata # generates a simple dataset the fits the expected API
sample_method(AnnData adata) -> AnnData adata # applies a simple modification that fits the method API
```

For adding datasets, methods and metrics, see above.

## Code Style and Testing

`singlecellopenproblems` is maintained at close to 100% code coverage. For datasets, methods, and metrics, tests are generated automatically. For additions outside this core functionality, contributors are encouraged to write tests for their code -- but if you do not know how to do so, please do not feel discouraged from contributing code! Others can always help you test your contribution.

Code style is dictated by [`black`](https://pypi.org/project/black/#installation-and-usage) and [`flake8`](https://flake8.pycqa.org/en/latest/) with [`hacking`](https://github.com/openstack/hacking). Code is automatically reformatted by [`pre-commit`](https://pre-commit.com/) when you push to GitHub.

## Code of Conduct

We abide by the principles of openness, respect, and consideration of others
of the Python Software Foundation: https://www.python.org/psf/codeofconduct/ and by
the Contributor Covenant. See our [Code of Conduct](CODE_OF_CONDUCT.md) for details.

## Attribution

This `CONTRIBUTING.md` was adapted from [scikit-learn](https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md).
