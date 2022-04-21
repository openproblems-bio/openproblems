from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import scanpy as sc


@method(
    method_name="Uniform Manifold Approximation and Projection (UMAP), "
    "as implemented by scanpy (logCPM, 1kHVG)",
    paper_name="UMAP: Uniform Manifold Approximation and Projection for "
    "Dimension Reduction",
    paper_url="https://arxiv.org/abs/1802.03426",
    paper_year=2018,
    code_url="https://github.com/lmcinnes/umap",
    code_version=check_version("umap-learn"),
)
def umap_logCPM_1kHVG(adata, test: bool = False, n_pca=50):
    adata = log_cpm_hvg(adata)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, use_rep="X_pca", n_pcs=n_pca)
    sc.tl.umap(adata)
    adata.obsm["X_emb"] = adata.obsm["X_umap"]
    return adata
