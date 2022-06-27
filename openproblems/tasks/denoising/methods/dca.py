from ....tools.decorators import method
from ....tools.utils import check_version

# import numpy as np
import scanpy as sc


def _dca(adata, test=False, epochs=None):
    if test:
        epochs = 30
    else:
        epochs = epochs or 300
    from dca.api import dca

    # find all-zero genes (columns)
    gene_sums = np.asarray(adata.obsm["train"].sum(axis=0)).flatten()
    is_missing = gene_sums == 0
    # make adata object with train counts
    adata2 = sc.AnnData(adata.obsm["train"])
    # mask all-zero genes
    adata2.X[:, is_missing] = 1
    # run DCA
    dca(adata2, epochs=epochs)
    adata.obsm["denoised"] = adata2.X  # adata2.X should call the count matrix of DCA.
    # return masked values to zero
    adata.obsm["denoised"][:.is_missing] = 0
    adata.uns["method_code_version"] = check_version("dca")
    return adata


@method(
    method_name="DCA",
    paper_name="Single-cell RNA-seq denoising using...",
    paper_url="https://www.nature.com/articles/s41467-018-07931-2",
    paper_year=2019,
    code_url="https://github.com/theislab/dca",
    image="openproblems-python-tf2.4",
)
def dca(adata, test=False, epochs=None):
    return _dca(adata, test=test, epochs=epochs)
