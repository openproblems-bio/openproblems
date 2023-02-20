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
    method_summary="TODO",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_reference="cable2021robust",
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
    # extract single cell reference data
    adata_sc, adata = split_sc_and_sp(adata)
    labels = np.unique(adata_sc.obs["label"])

    # set spatial coordinates for the single cell data
    adata_sc.obsm["spatial"] = np.ones((adata_sc.shape[0], 2))
    # remove rare cell types to prevent RCTD error
    celltype_counts = adata_sc.obs["label"].value_counts()
    adata_sc = adata_sc[
        ~adata_sc.obs["label"].isin(celltype_counts[celltype_counts < 25].index)
    ].copy()
    # run RCTD
    adata = _rctd(
        adata_sc, adata, fc_cutoff, fc_cutoff_reg, max_cores=multiprocessing.cpu_count()
    )

    # get predicted cell type proportions from obs
    cell_type_names = [f"xCT_{label}" for label in labels]

    # add proportions
    adata.obsm["proportions_pred"] = (
        adata.obs.reindex(cell_type_names, axis=1).fillna(0).to_numpy()
    )

    adata.uns["method_code_version"] = check_r_version("spacexr")

    return adata
