from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from ..utils import split_sc_and_sp
from typing import Optional

import pandas as pd

_seuratv3 = r_function("seuratv3.R", args="sce_sc, sce_sp, n_pcs")


@method(
    method_name="SeuratV3",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://doi.org/10.1016/j.cell.2019.05.031",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html",
    image="openproblems-r-extras",
)
def seuratv3(adata, test: bool = False, n_pca: Optional[int] = None):
    if test:
        n_pca = n_pca or 10
    else:  # pragma: nocover
        n_pca = n_pca or 30
    # extract single cell reference data
    adata_sc, adata = split_sc_and_sp(adata)
    # proportions_true gets lost in translation
    proportions_true = adata.obsm["proportions_true"]
    adata = _seuratv3(adata_sc, adata, n_pcs=n_pca)
    # get predicted cell type proportions from obs
    cell_type_names = pd.Index([x for x in adata.obs.columns if x.startswith("xCT_")])
    proportions_pred = adata.obs[cell_type_names].to_numpy()

    # add proportions
    adata.obsm["proportions_pred"] = proportions_pred
    adata.obsm["proportions_true"] = proportions_true

    adata.uns["method_code_version"] = check_r_version("Seurat")

    return adata
