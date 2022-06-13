import utils.warnings  # noqa: F401

# isort: split
import openproblems
import parameterized
import pytest
import utils.git
import utils.name

pytestmark = pytest.mark.skipif(
    len(utils.git.list_modified_tasks()) == 0, reason="No tasks have been modified"
)


@parameterized.parameterized.expand(
    [(metric,) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name.name_test,
)
def test_metric_metadata(metric):
    """Test for existence of metric metadata."""
    assert hasattr(metric, "metadata")
    for attr in ["metric_name", "maximize", "image"]:
        assert attr in metric.metadata
    assert isinstance(metric.metadata["maximize"], bool)


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            metric.__name__,
            metric.metadata["image"],
        )
        for task in utils.git.list_modified_tasks()
        for metric in task.METRICS
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test
def test_metric(task_name, metric_name, image):
    """Test computation of a metric."""
    import numbers

    task = getattr(openproblems.tasks, task_name)
    metric = getattr(task.metrics, metric_name)
    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
    )
    m = metric(adata)
    assert isinstance(m, numbers.Number)
