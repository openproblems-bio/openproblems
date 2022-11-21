from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import functools

_densmap_method = functools.partial(
    method,
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    image="openproblems-python-extras",
)


def _densmap(adata, obsm=None):
    from umap import UMAP

    if obsm:
        X = adata.obsm[obsm]
    else:
        X = adata.X
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(X)
    adata.uns["method_code_version"] = check_version("umap-learn")
    return adata


@_densmap_method(method_name="densMAP (logCPM, 1kHVG)")
def densmap_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    return _densmap(adata)


@_densmap_method(method_name="densMAP PCA (logCPM, 1kHVG)")
def densmap_pca_logCPM_1kHVG(adata, test: bool = False):
    import scanpy as sc

    adata = log_cpm_hvg(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    return _densmap(adata, obsm="X_pca")
