from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from .._utils import split_sc_and_sp

import pandas as pd

_seuratv3 = r_function("seuratv3.R")


@method(
    method_name="SeuratV3",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html",
    image="openproblems-r-extras",
)
def seuratv3(adata, test=False):
    # exctract single cell reference data
    sc_adata, sp_adata = split_sc_and_sp(adata)
    # store true proportions, due to error when passing DataFrame in obsm
    proportions_true = sp_adata.obsm["proportions_true"]
    # r function only accepts one argument, pass original adata object
    sp_adata = _seuratv3(adata)
    # get predicted cell type proportions from obs
    cell_type_names = pd.Index([x for x in sp_adata.obs.columns if "xCT_" in x[0:4]])
    proportions_pred = sp_adata.obs[cell_type_names]
    proportions_pred.columns = pd.Index([x.lstrip("xCT_") for x in cell_type_names])

    assert (
        proportions_pred.shape == proportions_true.shape
    ), "mismatched sizes between predicted and true proportions"

    # get observation order
    order = sp_adata.obs.index

    # add proportions
    sp_adata.obsm["proportions_pred"] = proportions_pred.iloc[order, :]
    sp_adata.obsm["proportions_true"] = proportions_true.iloc[order, :]

    sp_adata.uns["method_code_version"] = check_r_version("Seurat")

    return sp_adata
