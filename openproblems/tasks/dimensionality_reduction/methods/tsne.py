from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version


@method(
    method_name="“t-Distributed Stochastic Neighbor Embedding (t-SNE) (logCPM, 1kHVG)",
    paper_name="Visualizing Data using t-SNE",
    paper_url="https://www.jmlr.org/papers/v9/vandermaaten08a.html",
    paper_year=2008,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.manifold.TSNE.html#sklearn.manifold.TSNE",
    image="openproblems-python-extras",
)
def tsne_logCPM_1kHVG(adata, test: bool = False, n_pca=50):
    import scanpy as sc

    adata = log_cpm_hvg(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata, n_comps=n_pca, svd_solver="arpack")
    sc.tl.tsne(adata, use_rep="X_pca", n_pcs=n_pca)
    adata.obsm["X_emb"] = adata.obsm["X_tsne"]
    adata.uns["method_code_version"] = check_version("MulticoreTSNE")
    return adata
