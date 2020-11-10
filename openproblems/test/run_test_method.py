import sys
import anndata
import openproblems
from openproblems.test import utils

utils.ignore_warnings()


def test_method(task, dataset, method, data_path):
    adata = dataset(test=True)
    adata = method(adata)
    assert isinstance(adata, anndata.AnnData)
    assert task.checks.check_method(adata)
    adata.write_h5ad(data_path)


def main(task_name, method_name, dataset_name, data_path):
    openproblems.data.no_cleanup()
    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    method = getattr(task.methods, method_name)
    test_method(task, dataset, method, data_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
