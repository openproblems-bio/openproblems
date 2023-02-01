import openproblems
import os
import parameterized
import utils.asserts
import utils.docker
import utils.git
import utils.name

RETRIES = (
    int(os.environ["PYTEST_MAX_RETRIES"]) if "PYTEST_MAX_RETRIES" in os.environ else 2
)

METHOD_SUMMARY_MINLEN = 40
METHOD_SUMMARY_MAXLEN = 400


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            method.__name__,
            method.metadata["image"],
        )
        for task in openproblems.TASKS
        for method in task.METHODS
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test(timeout=600, retries=RETRIES)
def test_method(task_name, method_name, image):
    """Test application of a method."""
    import anndata
    import openproblems.utils

    task = getattr(openproblems.tasks, task_name)
    method = getattr(task.methods, method_name)
    adata = task.api.sample_dataset()
    openproblems.log.debug(
        "Testing {} method from {} task".format(method.__name__, task.__name__)
    )
    adata = method(adata, test=True)
    assert isinstance(adata, anndata.AnnData)
    assert task.api.check_method(adata, is_baseline=method.metadata["is_baseline"])
    if "method_code_version" not in adata.uns:
        openproblems.utils.future_warning(
            "Setting code_version in the method decorator is deprecated. "
            "Store code version in `adata.uns['method_code_version']` instead.",
            error_version="1.0",
            error_category=TypeError,
        )
        assert method.metadata["code_version"] is not None
    else:
        assert adata.uns["method_code_version"] != "ModuleNotFound"


@parameterized.parameterized.expand(
    [(method,) for task in openproblems.TASKS for method in task.METHODS],
    name_func=utils.name.name_test,
)
def test_method_metadata(method):
    """Test for existence of method metadata."""
    assert hasattr(method, "metadata")
    for attr in [
        "method_name",
        "method_summary",
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
    assert isinstance(method.metadata["method_summary"], str)
    assert len(method.metadata["method_summary"]) > METHOD_SUMMARY_MINLEN
    assert len(method.metadata["method_summary"]) < METHOD_SUMMARY_MAXLEN
    assert isinstance(method.metadata["paper_name"], str)
    assert isinstance(method.metadata["paper_year"], int)
    assert isinstance(method.metadata["paper_reference"], str)
    assert utils.asserts.assert_valid_reference(method.metadata["paper_reference"])
    assert isinstance(method.metadata["code_url"], str)
    assert utils.asserts.assert_url_accessible(method.metadata["code_url"])
    assert isinstance(method.metadata["is_baseline"], bool)
