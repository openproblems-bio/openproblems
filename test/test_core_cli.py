from openproblems.api.hash import docker_labels_from_api
from openproblems.api.hash import git_hash
from openproblems.api.main import main
from openproblems.api.utils import print_output

import importlib
import numpy as np
import openproblems
import os
import parameterized
import sklearn
import tempfile
import utils.name


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
    expected = np.array(
        [openproblems.utils.get_member_id(task) for task in openproblems.TASKS]
    )
    assert np.all(result == expected)
    result = np.array(main(["tasks"], do_print=True))
    expected = "\n".join(
        [openproblems.utils.get_member_id(task) for task in openproblems.TASKS] + [""]
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
            ["list", "--task", openproblems.utils.get_member_id(task), "--datasets"],
            do_print=False,
        )
    )
    expected = np.array([dataset.__name__ for dataset in task.DATASETS])
    assert np.all(result == expected)

    result = np.array(
        main(
            ["list", "--task", openproblems.utils.get_member_id(task), "--methods"],
            do_print=False,
        )
    )
    expected = np.array([method.__name__ for method in task.METHODS])
    assert np.all(result == expected)

    result = np.array(
        main(
            ["list", "--task", openproblems.utils.get_member_id(task), "--metrics"],
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
            openproblems.utils.get_member_id(task),
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


def test_hash_basic():
    assert main(["--test-hash"], do_print=False) is None


def test_version():
    assert main(["--version"], do_print=False) == openproblems.__version__


def test_help(capsys):
    assert main([], do_print=False) is None
    captured = capsys.readouterr()
    assert len(captured.out) > 0
    assert (
        "Open Problems for Single Cell Analysis command-line interface" in captured.out
    )


@parameterized.parameterized.expand(
    [
        ("label_projection", "--datasets", "pancreas_batch"),
        ("matching_modalities", "--methods", "mnn_log_scran_pooling"),
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


def test_hash_docker_api():
    assert docker_labels_from_api("circleci/python", tag="3.8-bullseye") is None
    labels = docker_labels_from_api("singlecellopenproblems/openproblems", tag="latest")
    assert "bio.openproblems.build" in labels
    assert "bio.openproblems.hash" in labels
    assert isinstance(labels["bio.openproblems.build"], str)
    assert isinstance(labels["bio.openproblems.hash"], str)
    assert labels["bio.openproblems.build"] in ["github_actions", "local"]


@parameterized.parameterized.expand(
    [
        (openproblems.tasks.label_projection.datasets.zebrafish_labs,),
        (openproblems.tasks.label_projection.methods.knn_classifier_log_cp10k,),
    ],
    name_func=utils.name.name_test,
)
def test_git_hash(func):
    h1 = git_hash(func)
    module = importlib.import_module(func.__wrapped__.__module__)
    assert git_hash(module) == h1
    assert git_hash(module.__file__) == h1


@parameterized.parameterized.expand(
    [
        (dataset, method, metric)
        for dataset in ["zebrafish_labs", None]
        for method in ["logistic_regression_log_cp10k", None]
        for metric in ["accuracy", None]
    ],
    name_func=utils.name.name_test,
)
def test_test(dataset, method, metric):
    """Test pipeline for dev testing."""
    args = ["test", "--task", "label_projection", "--test"]
    if dataset is not None:
        args += ["--dataset", dataset]
    if method is not None:
        args += ["--method", method]
    if metric is not None:
        args += ["--metric", metric]
    out = main(
        args,
        do_print=False,
    )
    if metric is None:
        assert out is None
    else:
        try:
            float(out)
        except ValueError:
            assert False, "result could not be converted to float"


def test_zero_metric():
    def __zero_metric(*args):
        return 0.0

    metric_name = utils.name.object_name(__zero_metric)
    task = openproblems.TASKS[0]
    setattr(task.metrics, metric_name, __zero_metric)
    adata = task.api.sample_dataset()
    with tempfile.TemporaryDirectory() as tempdir:
        dataset_file = os.path.join(tempdir, "dataset.h5ad")
        adata.var = adata.var[[]]
        adata.write_h5ad(dataset_file)
        result = main(
            [
                "evaluate",
                "--task",
                openproblems.utils.get_member_id(task),
                "--input",
                dataset_file,
                metric_name,
            ],
            do_print=True,
        )
        assert result == 0
        assert isinstance(result, int)


def test_pipeline():
    """Test evaluation pipeline."""
    with tempfile.TemporaryDirectory() as tempdir:
        dataset_file = os.path.join(tempdir, "dataset.h5ad")
        method_file = os.path.join(tempdir, "method.h5ad")
        version_file = os.path.join(tempdir, "version.txt")
        result_file = os.path.join(tempdir, "result.txt")
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
                "zebrafish_labs",
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
                "--version-file",
                version_file,
                "logistic_regression_log_cp10k",
            ],
            do_print=False,
        )
        assert os.path.isfile(method_file)
        assert os.path.isfile(version_file)
        with open(version_file, "r") as handle:
            code_version = handle.read()
        assert code_version == sklearn.__version__
        main(
            [
                "evaluate",
                "--task",
                "label_projection",
                "--input",
                method_file,
                "--output",
                result_file,
                "accuracy",
            ],
            do_print=False,
        )
        with open(result_file, "r") as handle:
            result = float(handle.read())
        assert isinstance(result, float)
