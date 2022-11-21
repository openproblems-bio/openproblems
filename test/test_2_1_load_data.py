import openproblems
import parameterized
import utils.docker
import utils.name


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            test,
            utils.TEMPDIR.name,
            dataset.metadata["image"],
        )
        for task in openproblems.TASKS
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
    utils.asserts.assert_array_equal(adata.X, adata.layers[task.DEFAULT_LAYER])
    adata2 = dataset(test=test)
    assert adata2.shape == adata.shape
    assert adata2.uns["_from_cache"]
    utils.cache.save(adata, tempdir, task, dataset, test=test)
