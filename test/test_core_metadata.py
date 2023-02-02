import openproblems
import openproblems.utils
import parameterized
import utils
import utils.asserts
import utils.cache
import utils.git
import utils.name

DATASET_SUMMARY_MINLEN = 40
DATASET_SUMMARY_MAXLEN = 1000


@parameterized.parameterized.expand(
    [(dataset,) for task in openproblems.TASKS for dataset in task.DATASETS],
    name_func=utils.name.name_test,
)
def test_dataset_metadata(dataset):
    """Test for existence of dataset metadata."""
    assert hasattr(dataset, "metadata")
    for attr in [
        "dataset_name",
        "data_url",
        "data_reference",
        "dataset_summary",
        "image",
    ]:
        assert attr in dataset.metadata
        assert dataset.metadata[attr] is not None

    assert isinstance(dataset.metadata["dataset_name"], str)
    assert isinstance(dataset.metadata["image"], str)
    assert dataset.metadata["image"].startswith("openproblems")
    assert isinstance(dataset.metadata["dataset_summary"], str)
    assert len(dataset.metadata["dataset_summary"]) > DATASET_SUMMARY_MINLEN
    assert len(dataset.metadata["dataset_summary"]) < DATASET_SUMMARY_MAXLEN
    assert isinstance(dataset.metadata["data_url"], str)
    assert utils.asserts.assert_url_accessible(dataset.metadata["data_url"])
    assert isinstance(dataset.metadata["data_reference"], str)
    assert utils.asserts.assert_valid_reference(dataset.metadata["data_reference"])


@parameterized.parameterized.expand(
    [(method,) for task in openproblems.TASKS for method in task.METHODS],
    name_func=utils.name.name_test,
)
def test_method_metadata(method):
    """Test for existence of method metadata."""
    assert hasattr(method, "metadata")
    for attr in [
        "method_name",
        "paper_name",
        "paper_reference",
        "paper_year",
        "code_url",
        "image",
        "is_baseline",
    ]:
        assert attr in method.metadata

    assert isinstance(method.metadata["image"], str)
    assert method.metadata["image"].startswith("openproblems")
    assert isinstance(method.metadata["method_name"], str)
    assert isinstance(method.metadata["paper_name"], str)
    assert isinstance(method.metadata["paper_year"], int)
    assert isinstance(method.metadata["paper_reference"], str)
    assert utils.asserts.assert_valid_reference(method.metadata["paper_reference"])
    assert isinstance(method.metadata["code_url"], str)
    assert utils.asserts.assert_url_accessible(method.metadata["code_url"])
    assert isinstance(method.metadata["is_baseline"], bool)


@parameterized.parameterized.expand(
    [(metric,) for task in openproblems.TASKS for metric in task.METRICS],
    name_func=utils.name.name_test,
)
def test_metric_metadata(metric):
    """Test for existence of metric metadata."""
    assert hasattr(metric, "metadata")
    for attr in ["metric_name", "maximize", "image"]:
        assert attr in metric.metadata
    assert isinstance(metric.metadata["maximize"], bool)
    assert isinstance(metric.metadata["metric_name"], str)
    assert isinstance(metric.metadata["image"], str)
    assert metric.metadata["image"].startswith("openproblems")
    assert isinstance(metric.metadata["paper_reference"], str)
    assert utils.asserts.assert_valid_reference(metric.metadata["paper_reference"])
