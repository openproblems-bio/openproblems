from . import utils

import sys
import openproblems
import anndata
import numbers

utils.ignore_warnings()


def test_metric(task, data_path, metric):
    """Test a metric applied to a specific dataset run through a specific method."""
    adata = anndata.read_h5ad(data_path)
    m = metric(adata)
    assert isinstance(m, numbers.Number)


def main(task_name, metric_name, data_path):
    """Run the test.

    This test is run via a subprocess call in test_metrics.py.
    """
    openproblems.data.no_cleanup()
    task = getattr(openproblems.tasks, task_name)
    metric = getattr(task.metrics, metric_name)
    openproblems.log.debug(
        "Testing {} metric on data located at {} from {} task".format(
            metric.__name__, data_path, task.__name__
        )
    )
    test_metric(task, data_path, metric)


if __name__ == "__main__":
    openproblems.log.debug("Running metric test with args {}".format(sys.argv))
    main(*sys.argv[1:])
