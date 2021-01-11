import utils
import parameterized
import openproblems

utils.warnings.ignore_warnings()


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            method.__name__,
            utils.TEMPDIR.name,
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
    import cache

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    adata = cache.load(
        tempdir, task, dataset, test=True, dependency="test_load_dataset"
    )
    openproblems.log.debug(
        "Testing {} method on {} dataset from {} task".format(
            method.__name__, dataset.__name__, task.__name__
        )
    )
    adata = method(adata)
    assert isinstance(adata, anndata.AnnData)
    assert task.checks.check_method(adata)
    cache.save(adata, tempdir, task, dataset, method=method)


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            method.__name__,
            metric.__name__,
            utils.TEMPDIR.name,
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
    """Test computation of a metric."""
    import cache
    import numbers

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    metric = getattr(task.metrics, metric_name)
    adata = cache.load(tempdir, task, dataset, method=method, dependency="test_method")
    openproblems.log.debug(
        "Testing {} metric on {} method applied to {} dataset from {} task".format(
            metric.__name__, method.__name__, dataset.__name__, task.__name__
        )
    )
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
