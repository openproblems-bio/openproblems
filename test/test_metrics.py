from . import utils

import parameterized
import openproblems

utils.ignore_warnings()


@parameterized.parameterized.expand(
    [(metric,) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name_test,
)
def test_metric_metadata(metric):
    """Test for existence of metric metadata."""
    assert hasattr(metric, "metadata")
    for attr in ["metric_name", "maximize", "image"]:
        assert attr in metric.metadata
    assert isinstance(metric.metadata["maximize"], bool)
