from ....tools.decorators import method
from ....tools.utils import check_version

# import numpy as np
import scanpy as sc


def _dca(adata):
    from dca.api import dca

    # sc.AnnData takes (counts, obs=[obs], vars=[vars]), but is tested to return an anndata even if just given counts. If DCA relies on obs or vars,
    # we will likely need to either access the vars of adata (vars index) and the barcodes (obs index) of
    # cells in adata.obs['train'], or else create fake matrices
    # by casting list(range(dim(adata.obsm['train'][1]) as header for vars and the next index [2] for header of obs
    adata2 = sc.AnnData(adata.obsm["train"])
    Y = dca(adata2, threads=1)
    adata.obsm["denoised"] = Y.X  # Y.X should call the count matrix of DCA.
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
def dca(adata, test=False):
    return _dca(adata)
