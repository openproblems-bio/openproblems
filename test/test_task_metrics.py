import openproblems
import parameterized
import utils.asserts
import utils.docker
import utils.name


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            metric.__name__,
            metric.metadata["image"],
        )
        for task in openproblems.TASKS
        for metric in task.METRICS
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test
def test_metric(task_name, metric_name, image):  # pragma: nocover
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
