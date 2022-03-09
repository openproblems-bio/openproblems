from ....tools.decorators import method


@method(
    method_name="KNN_Aggregate",
    image="openproblems-python-extras",
    paper_name="Chromatin Potential Identified by Shared Single-Cell Profiling\
                of RNA and Chromatin",
    paper_url="https://www.sciencedirect.com/science/article/pii/S0092867420312538",
    paper_year="2020",
    code_url=None,
    code_version="1.0",
)
def knn_aggregate(adata, n_neighbors=10):
    """Stabilize the ATAC profile by adding the signal of n_neighbors
    nearest neighbors for each cell.
    """

    import numpy as np
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="mode2_noisy", n_neighbors=n_neighbors)
    kernel = adata.obsp["distances"] > 0
    adata.obsm["mode2_denoised"] = np.dot(kernel, adata.obsm["mode2_noisy"])

    return adata
