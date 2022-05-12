import openproblems
import parameterized
import pytest
import utils.git
import utils.name
import utils.warnings

utils.warnings.ignore_warnings()
pytestmark = pytest.mark.skipif(
    len(utils.git.list_modified_tasks()) == 0, reason="No tasks have been modified"
)


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


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            metric.__name__,
            utils.TEMPDIR.name,
            metric.metadata["image"],
        )
        for task in utils.git.list_modified_tasks()
        for metric in task.METRICS
    ],
    name_func=utils.name.name_test,
    skip_on_empty=True,
)
@utils.docker.docker_test
def test_metric(task_name, metric_name, tempdir, image):
    """Test computation of a metric."""
    import numbers

    task = getattr(openproblems.tasks, task_name)
    metric = getattr(task.metrics, metric_name)
    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
    )
    m = metric(adata)
    assert isinstance(m, numbers.Number)


@parameterized.parameterized.expand(
    [
        (
            utils.TEMPDIR.name,
            openproblems.tasks.dimensionality_reduction.metrics.trustworthiness.metadata[  # noqa: E501
                "image"
            ],
        )
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_trustworthiness_sparse(tempdir, image):
    from scipy.sparse import csr_matrix

    task = openproblems.tasks.dimensionality_reduction
    metric = openproblems.tasks.dimensionality_reduction.metrics.trustworthiness

    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
    )
    adata.X = csr_matrix(adata.X)
    m = metric(adata)

    assert isinstance(m, float)
    assert 0 <= m <= 1


@parameterized.parameterized.expand(
    [
        (
            utils.TEMPDIR.name,
            openproblems.tasks.dimensionality_reduction.metrics.density_preservation.metadata[  # noqa: E501
                "image"
            ],
        )
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_density_preservation_matches_densmap(tempdir, image):
    from openproblems.tasks.dimensionality_reduction.metrics.density import _K
    from openproblems.tasks.dimensionality_reduction.metrics.density import _SEED
    from scipy.stats import pearsonr
    from umap import UMAP

    import numpy as np

    task = openproblems.tasks.dimensionality_reduction
    metric = openproblems.tasks.dimensionality_reduction.metrics.density_preservation

    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
    )

    (emb, ro, re) = UMAP(
        n_neighbors=_K, random_state=_SEED, densmap=True, output_dens=True
    ).fit_transform(adata.X)
    expected = pearsonr(ro, re)[0]

    adata.obsm["X_emb"] = emb
    actual = metric(adata)

    np.testing.assert_allclose(expected, actual)
