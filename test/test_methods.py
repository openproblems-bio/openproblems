import openproblems
import parameterized
import utils
import utils.docker
import utils.name
import utils.warnings

utils.warnings.ignore_warnings()


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            method.__name__,
            utils.TEMPDIR.name,
            method.metadata["image"],
        )
        for task in openproblems.TASKS
        for method in task.METHODS
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test(timeout=600, retries=2)
def test_method(task_name, method_name, tempdir, image):
    """Test application of a method."""
    import anndata
    import inspect

    task = getattr(openproblems.tasks, task_name)
    method = getattr(task.methods, method_name)
    adata = task.api.sample_dataset()
    openproblems.log.debug(
        "Testing {} method from {} task".format(method.__name__, task.__name__)
    )
    # use test mode of method for speed increase
    method_ = method
    if hasattr(method, "__wrapped__"):
        method_ = method.__wrapped__
    kwargs = {"test": True} if "test" in inspect.signature(method_).parameters else {}
    adata = method(adata, **kwargs)
    assert isinstance(adata, anndata.AnnData)
    assert task.api.check_method(adata)


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
        "code_version",
        "image",
    ]:
        assert attr in method.metadata
