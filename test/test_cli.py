from openproblems.api.main import main
from openproblems.api.utils import print_output

import numpy as np
import openproblems
import os
import parameterized
import tempfile
import utils


def test_print(capsys):
    """Test CLI output printer."""
    print_output("test")
    captured = capsys.readouterr()
    assert captured.out == "test\n"
    print_output(["foo", "bar"])
    captured = capsys.readouterr()
    assert captured.out == "foo\nbar\n"


def test_tasks(capsys):
    """Test task listing."""
    result = np.array(main(["tasks"], do_print=False))
    expected = np.array([task.__name__.split(".")[-1] for task in openproblems.TASKS])
    assert np.all(result == expected)
    result = np.array(main(["tasks"], do_print=True))
    expected = (
        "\n".join([task.__name__.split(".")[-1] for task in openproblems.TASKS]) + "\n"
    )
    captured = capsys.readouterr()
    assert captured.out == expected


@parameterized.parameterized.expand(
    [(task,) for task in openproblems.TASKS],
    name_func=utils.name.name_test,
)
def test_list(task):
    """Test function listing."""
    result = np.array(
        main(
            ["list", "--task", task.__name__.split(".")[-1], "--datasets"],
            do_print=False,
        )
    )
    expected = np.array([dataset.__name__ for dataset in task.DATASETS])
    assert np.all(result == expected)

    result = np.array(
        main(
            ["list", "--task", task.__name__.split(".")[-1], "--methods"],
            do_print=False,
        )
    )
    expected = np.array([method.__name__ for method in task.METHODS])
    assert np.all(result == expected)

    result = np.array(
        main(
            ["list", "--task", task.__name__.split(".")[-1], "--metrics"],
            do_print=False,
        )
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
        ],
        do_print=False,
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


@parameterized.parameterized.expand(
    [
        ("label_projection", "--datasets", "pancreas_batch"),
        ("multimodal_data_integration", "--methods", "mnn_log_scran_pooling"),
    ],
    name_func=utils.name.name_test,
)
def test_hash(task, function_type, function_name):
    """Test git hash function."""
    h1 = main(
        ["hash", "--task", task, function_type, function_name],
        do_print=False,
    )
    h2 = main(
        ["hash", "--task", task, function_type, function_name],
        do_print=False,
    )
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
                "pancreas_batch",
            ],
            do_print=False,
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
            ],
            do_print=False,
        )
        assert os.path.isfile(method_file)
        result = main(
            [
                "evaluate",
                "--task",
                "label_projection",
                "--input",
                method_file,
                "accuracy",
            ],
            do_print=False,
        )
        assert isinstance(result, float)
