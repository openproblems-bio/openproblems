import utils

import parameterized
import tempfile

import openproblems

utils.warnings.ignore_warnings()

TEMPDIR = tempfile.TemporaryDirectory()


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            method.__name__,
            TEMPDIR.name,
            method.metadata["image"],
        )
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for method in task.METHODS
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_method(task_name, dataset_name, method_name, tempdir, image):
    """Test application of a method."""
    import anndata
    import os

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    data_path = os.path.join(
        tempdir, "{}_{}_{}.h5ad".format(task_name, dataset_name, method_name)
    )
    openproblems.log.debug(
        "Testing {} method on {} dataset from {} task".format(
            method.__name__, dataset.__name__, task.__name__
        )
    )
    adata = dataset(test=True)
    adata = method(adata)
    assert isinstance(adata, anndata.AnnData)
    assert task.checks.check_method(adata)
    adata.write_h5ad(data_path)


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            method.__name__,
            metric.__name__,
            TEMPDIR.name,
            metric.metadata["image"],
        )
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for method in task.METHODS
        for metric in task.METRICS
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_metric(task_name, dataset_name, method_name, metric_name, tempdir, image):
    """Test application of a method."""
    import anndata
    import numbers
    import os

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    metric = getattr(task.metrics, metric_name)
    data_path = os.path.join(
        tempdir, "{}_{}_{}.h5ad".format(task_name, dataset_name, method_name)
    )
    openproblems.log.debug(
        "Testing {} metric on {} method applied to {} dataset from {} task".format(
            metric.__name__, method.__name__, dataset.__name__, task.__name__
        )
    )
    adata = anndata.read_h5ad(data_path)
    m = metric(adata)
    assert isinstance(m, numbers.Number)


@parameterized.parameterized.expand(
    [(method,) for task in openproblems.TASKS for method in task.METHODS],
    name_func=utils.name.name_test,
)
def test_method_metadata(method):
    """Test for existence of method metadata."""
    assert hasattr(method, "metadata")
    for attr in [
        "method_name",
        "paper_name",
        "paper_url",
        "paper_year",
        "code_url",
        "code_version",
        "image",
    ]:
        assert attr in method.metadata
