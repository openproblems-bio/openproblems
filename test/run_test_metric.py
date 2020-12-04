import utils

import sys
import openproblems
import anndata
import numbers
import nose2

utils.warnings.ignore_warnings()


def create_test(task, data_path, metric):
    """Test a metric applied to a specific dataset run through a specific method."""

    def test():
        adata = anndata.read_h5ad(data_path)
        m = metric(adata)
        assert isinstance(m, numbers.Number)

    return test


def main(task_name, metric_name, data_path):
    """Run the test.

    This test is run via a subprocess call in test_metrics.py.
    """
    global test_metric
    openproblems.data.no_cleanup()
    task = getattr(openproblems.tasks, task_name)
    metric = getattr(task.metrics, metric_name)
    openproblems.log.debug(
        "Testing {} metric on data located at {} from {} task".format(
            metric.__name__, data_path, task.__name__
        )
    )
    test_metric = create_test(task, data_path, metric)
    sys.argv = [
        "nose2",
        "--with-coverage",
        "--coverage-config",
        ".container_coveragerc",
    ]
    nose2.main()


if __name__ == "__main__":
    openproblems.log.debug("Running metric test with args {}".format(sys.argv))
    main(*sys.argv[1:])
