from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version

import functools

_seurat = r_function("seurat.R", args="sce, is_test")

_seurat_method = functools.partial(
    method,
    paper_name="Integrated analysis of multimodal single-cell data",
    paper_url="https://doi.org/10.1016/j.cell.2021.04.048",
    paper_year=2021,
    code_url="https://github.com/satijalab/seurat",
    image="openproblems-r-extras",
)


@_seurat_method(
    method_name="Seurat reference mapping (SCTransform)",
)
def seurat(adata, test=False):
    adata = _seurat(adata, test)
    adata.uns["method_code_version"] = check_r_version("Seurat")
    return adata
