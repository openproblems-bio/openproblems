from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import functools

_umap_method = functools.partial(
    method,
    paper_name="UMAP: Uniform Manifold Approximation and Projection for "
    "Dimension Reduction",
    paper_url="https://arxiv.org/abs/1802.03426",
    paper_year=2018,
    code_url="https://github.com/lmcinnes/umap",
)
_densmap_method = functools.partial(
    method,
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
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


@_umap_method(method_name="UMAP (logCPM, 1kHVG)")
def umap_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    return _umap(adata, genes=adata.var["highly_variable"])


@_umap_method(method_name="UMAP PCA (logCPM, 1kHVG)")
def umap_pca_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    return _umap(adata, n_comps=50, genes=adata.var["highly_variable"])


@_umap_method(method_name="UMAP (logCPM)")
def umap_logCPM(adata, test: bool = False):
    adata = log_cpm(adata)
    return _umap(adata)


@_umap_method(method_name="UMAP PCA (logCPM)")
def umap_pca_logCPM(adata, test: bool = False):
    adata = log_cpm(adata)
    return _umap(adata, n_comps=50)


@_densmap_method(method_name="densMAP (logCPM, 1kHVG)")
def densmap_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    return _umap(adata, densmap=True, genes=adata.var["highly_variable"])


@_densmap_method(method_name="densMAP PCA (logCPM, 1kHVG)")
def densmap_pca_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    return _umap(adata, densmap=True, n_comps=50, genes=adata.var["highly_variable"])


@_densmap_method(method_name="densMAP (logCPM)")
def densmap_logCPM(adata, test: bool = False):
    adata = log_cpm(adata)
    return _umap(adata, densmap=True)


@_densmap_method(method_name="densMAP PCA (logCPM)")
def densmap_pca_logCPM(adata, test: bool = False):
    adata = log_cpm(adata)
    return _umap(adata, densmap=True, n_comps=50)
