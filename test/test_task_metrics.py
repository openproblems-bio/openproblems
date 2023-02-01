import openproblems
import parameterized
import utils.asserts
import utils.docker
import utils.name

METRIC_SUMMARY_MINLEN = 40
METRIC_SUMMARY_MAXLEN = 400


@parameterized.parameterized.expand(
    [(metric,) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name.name_test,
)
def test_metric_metadata(metric):
    """Test for existence of metric metadata."""
    assert hasattr(metric, "metadata")
    for attr in ["metric_name", "metric_summary", "maximize", "image"]:
        assert attr in metric.metadata
    assert isinstance(metric.metadata["maximize"], bool)
    assert isinstance(metric.metadata["metric_name"], str)
    assert isinstance(metric.metadata["metric_summary"], str)
    assert len(metric.metadata["metric_summary"]) > METRIC_SUMMARY_MINLEN
    assert len(metric.metadata["metric_summary"]) < METRIC_SUMMARY_MAXLEN
    assert isinstance(metric.metadata["image"], str)
    assert metric.metadata["image"].startswith("openproblems")
    assert isinstance(metric.metadata["paper_reference"], str)
    assert utils.asserts.assert_valid_reference(metric.metadata["paper_reference"])


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
