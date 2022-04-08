from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import pandas as pd

_seuratv3 = r_function("seuratv3.R")


@method(
    method_name="SeuratV3",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html",
    code_version=check_version("seuratv3"),
    image="openproblems-r-extras",
)
def seuratv3(adata, test=False):
    # exctract single cell reference data
    sc_adata = adata.uns["sc_reference"].copy()
    # remove single cell reference from original anndata
    del adata.uns["sc_reference"]
    # set spatial coordinates for the single cell data
    sc_adata.obsm["spatial"] = np.ones((sc_adata.shape[0], 2))
    # add labels to spatial data to avoid NaN's in Seurat conversion
    adata.obs["label"] = -1 * np.ones(adata.shape[0])
    # store true proportions, due to error when passing DataFrame in obsm
    proportions_true = adata.obsm["proportions_true"]
    # concatenate single cell and spatial data, r_function only accepts one argument
    adata = adata.concatenate(
        sc_adata, batch_key="modality", batch_categories=["sp", "sc"]
    )
    # remove single cell reference anndata to reduce memory use
    del sc_adata
    # run SeuratV3
    adata = _seuratv3(adata)
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
