from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from typing import Optional

import functools
import pathlib

_seurat = r_function(
    "seurat_wrapper.R", args="sce, n_pcs, k_score=NULL, k_filter=NULL, script_path"
)

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
def seurat(
    adata,
    n_pcs: Optional[int] = None,
    k_score: Optional[int] = None,
    k_filter: Optional[int] = None,
    test: bool = False,
):
    kwargs = {}
    if test:
        kwargs["n_pcs"] = n_pcs or 5
        kwargs["k_score"] = k_score or 5
        kwargs["k_filter"] = k_filter or 20
    else:  # pragma: nocover
        kwargs["n_pcs"] = n_pcs or 50
        if k_score is not None:
            kwargs["k_score"] = k_score
        if k_filter is not None:
            kwargs["k_filter"] = k_filter
    adata = _seurat(
        adata,
        script_path=pathlib.Path(__file__).parent.joinpath("seurat.R").as_posix(),
        **kwargs,
    )
    adata.uns["method_code_version"] = check_r_version("Seurat")
    return adata
