from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.utils import check_version

import functools

_umap_method = functools.partial(
    method,
    paper_name=(
        "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction"
    ),
    paper_reference="mcinnes2018umap",
    paper_year=2018,
    code_url="https://github.com/lmcinnes/umap",
)
_densmap_method = functools.partial(
    method,
    paper_name=(
        "Assessing single-cell transcriptomic variability through"
        " density-preserving data visualization"
    ),
    paper_reference="narayan2021assessing",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    image="openproblems-python-extras",
)


def _umap(adata, n_comps=None, genes=None, densmap=False):
    from umap import UMAP

    import scanpy as sc

    if genes is not None:
        adata_input = adata[:, genes].copy()
    else:
        adata_input = adata

    if n_comps is not None:
        sc.tl.pca(adata_input, n_comps=n_comps, svd_solver="arpack")
        X = adata_input.obsm["X_pca"]
    else:
        X = adata_input.X

    adata.obsm["X_emb"] = UMAP(densmap=densmap, random_state=42).fit_transform(X)
    adata.uns["method_code_version"] = check_version("umap-learn")
    return adata


@_umap_method(method_name="UMAP (logCP10k, 1kHVG)")
def umap_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _umap(adata, genes=adata.var["highly_variable"])


@_umap_method(method_name="UMAP PCA (logCP10k, 1kHVG)")
def umap_pca_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _umap(adata, n_comps=50, genes=adata.var["highly_variable"])


@_umap_method(method_name="UMAP (logCP10k)")
def umap_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _umap(adata)


@_umap_method(method_name="UMAP PCA (logCP10k)")
def umap_pca_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _umap(adata, n_comps=50)


@_densmap_method(method_name="densMAP (logCP10k, 1kHVG)")
def densmap_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _umap(adata, densmap=True, genes=adata.var["highly_variable"])


@_densmap_method(method_name="densMAP PCA (logCP10k, 1kHVG)")
def densmap_pca_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _umap(adata, densmap=True, n_comps=50, genes=adata.var["highly_variable"])


@_densmap_method(method_name="densMAP (logCP10k)")
def densmap_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _umap(adata, densmap=True)


@_densmap_method(method_name="densMAP PCA (logCP10k)")
def densmap_pca_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _umap(adata, densmap=True, n_comps=50)
