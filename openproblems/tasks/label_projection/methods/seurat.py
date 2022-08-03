from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version

import functools

_seurat = r_function("seurat.R")

_seurat_method = functools.partial(
    method,
    paper_name="Integrated analysis of multimodal single-cell data",
    paper_url="https://www.sciencedirect.com/science/article/pii/S0092867421005833",
    paper_year=2021,
    code_url="https://github.com/satijalab/seurat",
    image="openproblems-r-extras",
)

@_seurat_method(
    method_name="Seurat reference mapping (SCTransform)",
)
def seurat(adata, test=False):
    adata = _seurat(adata)
    adata.uns["method_code_version"] = check_r_version("Seurat")
    return adata