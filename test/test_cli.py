from . import utils
import numpy as np
import openproblems
from openproblems.api.main import main
import os
import parameterized
import tempfile


def test_tasks():
    """Test task listing."""
    result = np.array(main(["tasks"]))
    expected = np.array([task.__name__.split(".")[-1] for task in openproblems.TASKS])
    assert np.all(result == expected)


@parameterized.parameterized.expand(
    [(task,) for task in openproblems.TASKS],
    name_func=utils.name.name_test,
)
def test_list(task):
    """Test function listing."""
    result = np.array(
        main(["list", "--task", task.__name__.split(".")[-1], "--datasets"])
    )
    expected = np.array([dataset.__name__ for dataset in task.DATASETS])
    assert np.all(result == expected)

    result = np.array(
        main(["list", "--task", task.__name__.split(".")[-1], "--methods"])
    )
    expected = np.array([method.__name__ for method in task.METHODS])
    assert np.all(result == expected)

    result = np.array(
        main(["list", "--task", task.__name__.split(".")[-1], "--metrics"])
    )
    expected = np.array([metric.__name__ for metric in task.METRICS])
    assert np.all(result == expected)


def _test_image(task, function_type, function):
    result = main(
        [
            "image",
            "--task",
            task.__name__.split(".")[-1],
            function_type,
            function.__name__,
        ]
    )
    expected = function.metadata["image"]
    assert result == expected


@parameterized.parameterized.expand(
    [(task, dataset) for task in openproblems.TASKS for dataset in task.DATASETS],
    name_func=utils.name.name_test,
)
def test_image_datasets(task, dataset):
    """Test Docker image retrieval for datasets."""
    _test_image(task, "--datasets", dataset)


@parameterized.parameterized.expand(
    [(task, method) for task in openproblems.TASKS for method in task.METHODS],
    name_func=utils.name.name_test,
)
def test_image_methods(task, method):
    """Test Docker image retrieval for methods."""
    _test_image(task, "--methods", method)


@parameterized.parameterized.expand(
    [(task, metric) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name.name_test,
)
def test_image_metrics(task, metric):
    """Test Docker image retrieval for metrics."""
    _test_image(task, "--metrics", metric)


def test_hash():
    """Test git hash function."""
    h1 = main(["hash", "--task", "label_projection", "--datasets", "pancreas"])
    h2 = main(["hash", "--task", "label_projection", "--datasets", "pancreas"])
    assert h1 == h2


def test_pipeline():
    """Test evaluation pipeline."""
    with tempfile.TemporaryDirectory() as tempdir:
        dataset_file = os.path.join(tempdir, "dataset.h5ad")
        method_file = os.path.join(tempdir, "method.h5ad")
        assert not os.path.isfile(dataset_file)
        assert not os.path.isfile(method_file)
        main(
            [
                "load",
                "--task",
                "label_projection",
                "--test",
                "--output",
                dataset_file,
                "pancreas",
            ]
        )
        assert os.path.isfile(dataset_file)
        main(
            [
                "run",
                "--task",
                "label_projection",
                "--input",
                dataset_file,
                "--output",
                method_file,
                "logistic_regression_log_cpm",
            ]
        )
        assert os.path.isfile(method_file)
        result = main(
            [
                "evaluate",
                "--task",
                "label_projection",
                "--input",
                dataset_file,
                "--output",
                method_file,
                "accuracy",
            ]
        )
        assert isinstance(result, float)
