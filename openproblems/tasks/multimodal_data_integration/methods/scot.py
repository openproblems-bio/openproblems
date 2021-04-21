from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version

import sklearn.decomposition


def _scot(adata, n_svd=100, balanced=False):
    from SCOT import SCOT

    # PCA reduction
    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])

    # Initialize SCOT
    scot = SCOT(X_pca, Y_pca)

    # call the unbalanced alignment
    # From https://github.com/rsinghlab/SCOT/blob/master/examples/unbalanced_GW_SNAREseq.ipynb # noqa: 501
    X_new_unbal, y_new_unbal = scot.align(
        k=50, e=1e-3, rho=0.0005, normalize=True, balanced=balanced
    )
    adata.obsm["aligned"] = X_new_unbal
    adata.obsm["mode2_aligned"] = y_new_unbal

    return adata


@method(
    method_name="Single Cell Optimal Transport (sqrt CPM unbalanced)",
    paper_name="Gromov-Wasserstein optimal transport"
    "to align single-cell multi-omics data",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.04.28.066787",
    paper_year=2020,
    code_url="https://github.com/rsinghlab/SCOT",
    code_version=check_version("SCOT"),
    image="openproblems-python-extras",
)
def scot_sqrt_cpm_unbalanced(adata, n_svd=100, balanced=False):
    sqrt_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _scot(adata, n_svd=n_svd, balanced=balanced)
    return adata


@method(
    method_name="Single Cell Optimal Transport (sqrt CPM balanced)",
    paper_name="Gromov-Wasserstein optimal transport to "
    "align single-cell multi-omics data",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.04.28.066787",
    paper_year=2020,
    code_url="https://github.com/rsinghlab/SCOT",
    code_version=check_version("SCOT"),
    image="openproblems-python-extras",
)
def scot_sqrt_cpm_balanced(adata, n_svd=100, balanced=True):
    sqrt_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _scot(adata, n_svd=n_svd, balanced=balanced)
    return adata


@method(
    method_name="Single Cell Optimal Transport (log scran unbalanced)",
    paper_name="Gromov-Wasserstein optimal transport to "
    "align single-cell multi-omics data",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.04.28.066787",
    paper_year=2020,
    code_url="https://github.com/rsinghlab/SCOT",
    code_version=check_version("SCOT"),
    image="openproblems-python-extras",
)
def scot_log_scran_pooling_unbalanced(adata, n_svd=100, balanced=False):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _scot(adata, n_svd=n_svd, balanced=balanced)
    return adata


@method(
    method_name="Single Cell Optimal Transport (log scran balanced)",
    paper_name="Gromov-Wasserstein optimal transport to "
    "align single-cell multi-omics data",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.04.28.066787",
    paper_year=2020,
    code_url="https://github.com/rsinghlab/SCOT",
    code_version=check_version("SCOT"),
    image="openproblems-python-extras",
)
def scot_log_scran_pooling_balanced(adata, n_svd=100, balanced=True):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _scot(adata, n_svd=n_svd, balanced=balanced)
    return adata
