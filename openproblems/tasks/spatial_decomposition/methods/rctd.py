from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from ..utils import split_sc_and_sp

import numpy as np

_rctd = r_function("rctd.R", args="sce_sc, sce_sp")

RCTD_MIN_CELLTYPE_COUNT = 25


@method(
    method_name="RCTD",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_url="https://doi.org/10.1038/s41587-021-00830-w",
    paper_year=2020,
    code_url="https://github.com/dmcable/spacexr",
    image="openproblems-r-extras",
)
def rctd(adata, test=False):
    # exctract single cell reference data
    adata_sc, adata = split_sc_and_sp(adata)

    # set spatial coordinates for the single cell data
    adata_sc.obsm["spatial"] = np.ones((adata_sc.shape[0], 2))
    # run RCTD
    adata = _rctd(adata_sc, adata)

    # get predicted cell type proportions from obs
    cell_type_names = [x for x in adata.obs.columns if x.startswith("xCT")]

    # add proportions
    adata.obsm["proportions_pred"] = adata.obs[cell_type_names].to_numpy()

    adata.uns["method_code_version"] = check_r_version("spacexr")

    return adata
