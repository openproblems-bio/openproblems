import sys
import openproblems


def test_method(task, dataset, method, data_path):
    adata = dataset(test=True)
    output = method(adata)
    assert output is None
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
