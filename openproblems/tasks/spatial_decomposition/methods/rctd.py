from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np

_rctd = r_function("rctd.R")


@method(
    method_name="RCTD",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_url="https://www.nature.com/articles/s41587-021-00830-w",
    paper_year=2020,
    code_url="https://github.com/almaan/RCTD",
    code_version=check_version("rctd"),
)
def rctd(adata):
    # exctract single cell reference data
    sc_adata = adata.uns["sc_reference"].copy()
    # remove single cell reference from original anndata
    del adata.uns["sc_reference"]
    # set spatial coordinates for the single cell data
    sc_adata.obsm["spatial"] = np.ones((sc_adata.shape[0], 2))
    # concatenate single cell and spatial data, r_function only accepts one argument
    adata = adata.concatenate(
        sc_adata, batch_key="modality", batch_categories=["sp", "sc"]
    )
    # remove single cell reference anndata to reduce memory use
    del sc_adata
    # run RCTD
    adata = _rctd(adata)
    return adata
