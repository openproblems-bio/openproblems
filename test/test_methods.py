import openproblems
import parameterized
import pytest
import utils.docker
import utils.git
import utils.name
import utils.warnings

utils.warnings.ignore_warnings()
pytestmark = pytest.mark.skipif(
    len(utils.git.list_modified_tasks()) == 0, reason="No tasks have been modified"
)


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            method.__name__,
            utils.TEMPDIR.name,
            method.metadata["image"],
        )
        for task in utils.git.list_modified_tasks()
        for method in task.METHODS
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test(timeout=600, retries=2)
def test_method(task_name, method_name, tempdir, image):
    """Test application of a method."""
    import anndata

    task = getattr(openproblems.tasks, task_name)
    method = getattr(task.methods, method_name)
    adata = task.api.sample_dataset()
    openproblems.log.debug(
        "Testing {} method from {} task".format(method.__name__, task.__name__)
    )
    adata = method(adata, test=True)
    assert isinstance(adata, anndata.AnnData)
    assert task.api.check_method(adata)
    if "method_code_version" not in adata.uns:
        utils.warnings.future_warning(
            "Setting code_version in the method decorator is deprecated. "
            "Store code version in `adata.uns['method_code_version']` instead.",
            error_version="1.0",
            error_category=TypeError,
        )
        assert method.metadata["code_version"] is not None


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
        "paper_url",
        "paper_year",
        "code_url",
        "image",
    ]:
        assert attr in method.metadata
