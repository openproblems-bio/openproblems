from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from ..utils import split_sc_and_sp
from typing import Optional

import multiprocessing
import numpy as np

_rctd = r_function("rctd.R", args="sce_sc, sce_sp, fc_cutoff, fc_cutoff_reg, max_cores")


@method(
    method_name="RCTD",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_url="https://doi.org/10.1038/s41587-021-00830-w",
    paper_year=2020,
    code_url="https://github.com/dmcable/spacexr",
    image="openproblems-r-extras",
)
def rctd(
    adata,
    fc_cutoff: Optional[float] = None,
    fc_cutoff_reg: Optional[float] = None,
    test=False,
):
    if test:
        fc_cutoff = fc_cutoff or 0.05
        fc_cutoff_reg = fc_cutoff_reg or 0.075
    else:  # pragma: nocover
        fc_cutoff = fc_cutoff or 0.5
        fc_cutoff_reg = fc_cutoff_reg or 0.75
    # exctract single cell reference data
    adata_sc, adata = split_sc_and_sp(adata)

    # set spatial coordinates for the single cell data
    adata_sc.obsm["spatial"] = np.ones((adata_sc.shape[0], 2))
    # run RCTD
    adata = _rctd(
        adata_sc, adata, fc_cutoff, fc_cutoff_reg, max_cores=multiprocessing.cpu_count()
    )

    # get predicted cell type proportions from obs
    cell_type_names = [x for x in adata.obs.columns if x.startswith("xCT")]

    # add proportions
    adata.obsm["proportions_pred"] = adata.obs[cell_type_names].to_numpy()

    adata.uns["method_code_version"] = check_r_version("spacexr")

    return adata
