"""Specific tests for the dimensionality_reduction task"""
import numpy as np
import openproblems
import parameterized
import utils.docker
import utils.git

# global skip
TASK = openproblems.tasks.dimensionality_reduction


@utils.docker.docker_test(image=TASK.metrics.trustworthiness.metadata["image"])
def test_trustworthiness_sparse():  # pragma: nocover
    from scipy.sparse import csr_matrix

    metric = TASK.metrics.trustworthiness

    adata = TASK.api.sample_dataset()
    adata = TASK.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, TASK.__name__)
    )
    adata.X = csr_matrix(adata.X)
    m = metric(adata)

    assert isinstance(m, float)
    assert 0 <= m <= 1


def test_density_preservation_matches_densmap():
    from openproblems.tasks.dimensionality_reduction.metrics.density import _K
    from openproblems.tasks.dimensionality_reduction.metrics.density import _SEED
    from scipy.stats import pearsonr
    from umap import UMAP

    import numpy as np

    metric = TASK.metrics.density_preservation

    adata = TASK.api.sample_dataset()
    adata = TASK.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, TASK.__name__)
    )

    (emb, ro, re) = UMAP(
        n_neighbors=_K, random_state=_SEED, densmap=True, output_dens=True
    ).fit_transform(adata.X)
    expected = pearsonr(ro, re)[0]

    adata.obsm["X_emb"] = emb
    actual = metric(adata)

    np.testing.assert_allclose(expected, actual, rtol=1e-3)


@parameterized.parameterized.expand(
    [(200,), (1000,)],
    name_func=utils.name.name_test,
)
def test_distance_correlation_with_svd(n_svd):
    import numpy as np

    metric = TASK.metrics.distance_correlation

    adata = TASK.api.sample_dataset()
    adata = TASK.api.sample_method(adata)
    adata.obsm["X_emb"] = adata.X.toarray()

    expected = 1
    actual = metric(adata, n_svd=n_svd)

    np.testing.assert_allclose(expected, actual, rtol=1e-3)


def test_density_preservation_perfect():
    import numpy as np

    metric = TASK.metrics.density_preservation

    adata = TASK.api.sample_dataset()
    adata = TASK.api.sample_method(adata)

    adata.obsm["X_emb"] = adata.X.toarray()
    actual = metric(adata)

    np.testing.assert_allclose(1, actual)


def test_diffusion_map_no_convergence():
    import numpy as np
    import scipy.sparse.linalg

    adata = TASK.datasets.olsson_2016_mouse_blood()
    # no exception with retries
    adata = TASK.methods.diffusion_map(adata)
    # exception with no retries
    np.testing.assert_raises(
        scipy.sparse.linalg.ArpackNoConvergence,
        openproblems.tasks.dimensionality_reduction.methods.diffusion_map,
        adata,
        n_retries=0,
    )


@parameterized.parameterized.expand(
    [("qlocal",), ("qglobal",)],
    name_func=utils.name.name_test,
)
def test_pydrmetrics_subsample(metric_name):
    np.random.seed(0)
    metric = getattr(openproblems.tasks.dimensionality_reduction.metrics, metric_name)

    adata = TASK.api.sample_dataset()
    adata = TASK.api.sample_method(adata)
    metric_full = metric(adata)
    metric_sub = metric(adata, max_samples=int(adata.shape[0] * 0.9))
    np.testing.assert_allclose(metric_full, metric_sub, atol=0.05)
