from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import scanpy as sc


@method(
    method_name="“t-Distributed Stochastic Neighbor Embedding (t-SNE) (logCPM, 1kHVG)",
    paper_name="Visualizing Data using t-SNE",
    paper_url="https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf",
    paper_year=2008,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.manifold.TSNE.html#sklearn.manifold.TSNE",
    code_version=check_version("MulticoreTSNE"),
    image="openproblems-python-extras",
)
def tsne_logCPM_1kHVG(adata, test: bool = False, n_pca=50):
    adata = log_cpm_hvg(adata)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    sc.tl.tsne(adata, use_rep="X_pca", n_pcs=n_pca)
    adata.obsm["X_emb"] = adata.obsm["X_tsne"]
    return adata
