from ....tools.decorators import method
from ....tools.utils import check_version
from .preprocessing import preprocess_logCPM_1kHVG

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
def umap_logCPM_1kHVG(adata):
    preprocess_logCPM_1kHVG(adata)
    sc.pp.neighbors(adata, use_rep="X_input", n_pcs=50)
    sc.tl.umap(adata)
    adata.obsm["X_emb"] = adata.obsm["X_umap"]
    return adata
