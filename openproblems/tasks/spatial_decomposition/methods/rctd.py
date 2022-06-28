from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import split_sc_and_sp

import numpy as np

_rctd = r_function("rctd.R", args="sce_sc, sce_sp")


@method(
    method_name="RCTD",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_url="https://www.nature.com/articles/s41587-021-00830-w",
    paper_year=2020,
    code_url="https://github.com/almaan/RCTD",
    code_version=check_version("rctd"),
    image="openproblems-r-extras",
)
def rctd(adata, test=False):
    # exctract single cell reference data
    adata_sc, adata_sp = split_sc_and_sp(adata)
    del adata
    # set spatial coordinates for the single cell data
    adata_sc.obsm["spatial"] = np.ones((adata_sc.shape[0], 2))
    # run RCTD
    adata_sp = _rctd(adata_sc, adata_sp)

    # get predicted cell type proportions from obs
    cell_type_names = [x for x in adata_sp.obs.columns if x.startswith("xCT")]

    # add proportions
    adata_sp.obsm["proportions_pred"] = adata_sp.obs[cell_type_names].to_numpy()

    return adata_sp
