import parameterized

import openproblems
from openproblems.test import utils

utils.ignore_warnings()


@parameterized.parameterized.expand(
    [(metric,) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name_test,
)
def test_metric_metadata(metric):
    assert hasattr(metric, "metadata")
    for attr in ["metric_name", "maximize", "image"]:
        assert attr in metric.metadata
    assert isinstance(metric.metadata["maximize"], bool)
