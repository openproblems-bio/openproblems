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
            dataset.__name__,
            test,
            utils.TEMPDIR.name,
            dataset.metadata["image"],
        )
        for task in utils.git.list_modified_tasks()
        for dataset in task.DATASETS
        for test in [True]
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test(retries=2)
def test_load_dataset(task_name, dataset_name, test, tempdir, image):
    """Test loading and caching of a dataset."""
    import utils.asserts
    import utils.cache

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    adata = dataset(test=test)
    utils.asserts.assert_finite(adata.X)
    utils.asserts.assert_array_equal(adata.X, adata.layers["counts"])
    adata2 = dataset(test=test)
    assert adata2.shape == adata.shape
    assert adata2.uns["_from_cache"]
    utils.cache.save(adata, tempdir, task, dataset, test=test)
