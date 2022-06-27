from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import split_sc_and_sp

import numpy as np
import pandas as pd

_rctd = r_function("rctd.R")


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
    sc_adata, adata = split_sc_and_sp(adata)
    # set spatial coordinates for the single cell data
    sc_adata.obsm["spatial"] = np.ones((sc_adata.shape[0], 2))
    # store true proportions, due to error when passing DataFrame in obsm
    proportions_true = adata.obsm["proportions_true"]
    # concatenate single cell and spatial data, r_function only accepts one argument
    adata = adata.concatenate(
        sc_adata, batch_key="modality", batch_categories=["sp", "sc"]
    )
    # remove single cell reference anndata to reduce memory use
    del sc_adata
    # run RCTD
    adata = _rctd(adata)
    # remove appended subsetting name from index names
    new_idx = pd.Index([x.rstrip("-sp") for x in adata.obs.index], name="sample")
    adata.obs_names = new_idx
    adata.obsm.dim_names = new_idx

    # get predicted cell type proportions from obs
    cell_type_names = [x for x in adata.obs.columns if "xCT_" in x[0:4]]
    proportions_pred = adata.obs[cell_type_names]
    proportions_pred.columns = pd.Index([x.lstrip("xCT_") for x in cell_type_names])
    # make sure true proportions are concurrently ordered with predicted
    proportions_true = proportions_true.loc[new_idx, :]

    # add proportions
    adata.obsm["proportions_pred"] = proportions_pred
    adata.obsm["proportions_true"] = proportions_true

    return adata
