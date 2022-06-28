from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from .._utils import split_sc_and_sp

import pandas as pd

_seuratv3 = r_function("seuratv3.R", args="sce_sc, sce_sp")


@method(
    method_name="SeuratV3",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html",
    image="openproblems-r-extras",
)
def seuratv3(adata, test=False):
    # extract single cell reference data
    adata_sc, adata = split_sc_and_sp(adata)
    # proportions_true gets lost in translation
    proportions_true = adata.obsm["proportions_true"]
    adata = _seuratv3(adata_sc, adata)
    # get predicted cell type proportions from obs
    cell_type_names = pd.Index([x for x in adata.obs.columns if x.startswith("xCT_")])
    proportions_pred = adata.obs[cell_type_names].to_numpy()

    # add proportions
    adata.obsm["proportions_pred"] = proportions_pred
    adata.obsm["proportions_true"] = proportions_true

    adata.uns["method_code_version"] = check_r_version("Seurat")

    return adata
