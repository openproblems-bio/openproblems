from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import functools

_tsne_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name="Visualizing Data using t-SNE",
    paper_reference="vandermaaten2008visualizing",
    paper_year=2008,
    code_url=(
        "https://scikit-learn.org/stable/modules/generated/"
        "sklearn.manifold.TSNE.html#sklearn.manifold.TSNE"
    ),
    image="openproblems-python-extras",
)


def _tsne(adata, genes=None, test=False, n_pca=50):
    import scanpy as sc

    if genes is not None:
        X = adata[:, genes].copy().X
    else:
        X = adata.X

    adata.obsm["X_pca"] = sc.tl.pca(X, n_comps=n_pca, svd_solver="arpack")
    sc.tl.tsne(adata, use_rep="X_pca", n_pcs=n_pca)
    adata.obsm["X_emb"] = adata.obsm["X_tsne"]
    adata.uns["method_code_version"] = check_version("MulticoreTSNE")
    return adata


@_tsne_method(
    method_name="t-Distributed Stochastic Neighbor Embedding (t-SNE) (logCPM, 1kHVG)"
)
def tsne_logCPM_1kHVG(adata, test: bool = False, n_pca=50):
    adata = log_cpm_hvg(adata)
    return _tsne(adata, genes=adata.var["highly_variable"], test=test, n_pca=n_pca)


@_tsne_method(
    method_name="t-Distributed Stochastic Neighbor Embedding (t-SNE) (logCPM)"
)
def tsne_logCPM(adata, test: bool = False, n_pca=50):
    adata = log_cpm(adata)
    return _tsne(adata, test=test, n_pca=n_pca)
