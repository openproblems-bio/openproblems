import sys
import anndata
import openproblems
from openproblems.test import utils, log_level

utils.ignore_warnings()


def test_method(task, dataset, method, data_path):
    """Test a method applied to a specific dataset."""
    adata = dataset(test=True)
    adata = method(adata)
    assert isinstance(adata, anndata.AnnData)
    assert task.checks.check_method(adata)
    adata.write_h5ad(data_path)


def main(task_name, method_name, dataset_name, data_path):
    """Run the test.

    This test is run via a subprocess call in test_methods.py.
    """
    openproblems.data.no_cleanup()
    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    openproblems.log.debug(
        "Testing {} method on {} dataset from {} task".format(
            method.__name__, dataset.__name__, task.__name__
        )
    )
    test_method(task, dataset, method, data_path)


if __name__ == "__main__":
    openproblems.log.setLevel(log_level)
    openproblems.log.debug("Running method test with args {}".format(sys.argv))
    main(*sys.argv[1:])
