"""Specific tests for the dimensionality_reduction task"""
import openproblems
import utils.docker
import utils.git

# global skip
TASK = openproblems.tasks.dimensionality_reduction


@utils.docker.docker_test(image=TASK.metrics.trustworthiness.metadata["image"])
def test_trustworthiness_sparse():  # pragma: nocover
    from scipy.sparse import csr_matrix

    task = openproblems.tasks.dimensionality_reduction
    metric = task.metrics.trustworthiness

    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)
    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
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

    np.testing.assert_allclose(expected, actual, rtol=1e-3)


def test_density_preservation_perfect():
    import numpy as np

    task = openproblems.tasks.dimensionality_reduction
    metric = openproblems.tasks.dimensionality_reduction.metrics.density_preservation

    adata = task.api.sample_dataset()
    adata = task.api.sample_method(adata)

    adata.obsm["X_emb"] = adata.X.toarray()
    actual = metric(adata)

    np.testing.assert_allclose(1, actual)


def test_diffusion_map_no_convergence():
    import numpy as np
    import scipy.sparse.linalg

    adata = (
        openproblems.tasks.dimensionality_reduction.datasets.olsson_2016_mouse_blood()
    )
    # no exception with retries
    adata = openproblems.tasks.dimensionality_reduction.methods.diffusion_map(adata)
    # exception with no retries
    np.testing.assert_raises(
        scipy.sparse.linalg.ArpackNoConvergence,
        openproblems.tasks.dimensionality_reduction.methods.diffusion_map,
        adata,
        n_retries=0,
    )
