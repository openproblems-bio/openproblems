from ....tools.decorators import method

import scanpy as sc

@method(
    method_name="Uniform Manifold Approximation and Projection (UMAP)",
    paper_name="UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction",
    paper_url="https://arxiv.org/abs/1802.03426",
    paper_year=2018,
    code_url="https://github.com/lmcinnes/umap",
    code_version="0.4.6",
)


def umap(adata):
    
    sc.tl.umap(adata)
