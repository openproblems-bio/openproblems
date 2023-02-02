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
def test_method(task_name, method_name, image):  # pragma: nocover
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
            "Setting code_version in the method decorator is deprecated. Store code"
            " version in `adata.uns['method_code_version']` instead.",
            error_version="1.0",
            error_category=TypeError,
        )
        assert method.metadata["code_version"] is not None
    else:
        assert adata.uns["method_code_version"] != "ModuleNotFound"
