# SingleCellOpenProblems

[![Travis CI Build](https://api.travis-ci.com/scottgigante/SingleCellOpenProblems.svg?branch=master)](https://travis-ci.com/scottgigante/SingleCellOpenProblems)
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Formalizing and benchmarking open problems in single-cell genomics

[Results](results.md)

## API

Each task consists of datasets, methods, and metrics.

Datasets should take no arguments and return an AnnData object.

```
function dataset() -> AnnData adata
```

Methods should take an AnnData object and store the output in-place in  `adata.obs` according to the specification of the task.

```
function method(AnnData adata) -> None
```

Metrics should take an AnnData object and return a float.

```
function metric(AnnData adata) -> float
```

Task-specific APIs are described in the README for each task.

* [Label Projection](openproblems/tasks/label_projection)

## Adding a new dataset / method / metric

To add a new dataset, method, or metric to a task, simply create a new `.py` file corresponding to your proposed new functionality and import the main function in the corresponding `__init__.py`. E.g., to add a "F2" metric to the label projection task, we would create `openproblems/tasks/label_projection/metrics/f2.py` and add a line 
```
from .f2 import f2
```
to [`openproblems/tasks/label_projection/metrics/__init__.py`](openproblems/tasks/label_projection/metrics/__init__.py).
