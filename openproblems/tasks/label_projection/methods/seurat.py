from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_r_version
from typing import Optional

import pathlib

_seurat = r_function(
    "seurat_wrapper.R", args="sce, n_pcs, k_score, k_filter, script_path"
)


@method(
    method_name="Seurat reference mapping (SCTransform)",
    paper_name="Integrated analysis of multimodal single-cell data",
    paper_reference="hao2021integrated",
    paper_year=2021,
    code_url="https://github.com/satijalab/seurat",
    image="openproblems-r-extras",
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
        kwargs["k_score"] = k_score or 30
        kwargs["k_filter"] = k_filter or 200
    adata = _seurat(
        adata,
        script_path=pathlib.Path(__file__).parent.joinpath("seurat.R").as_posix(),
        **kwargs,
    )
    # R conversion scrambles booleans
    adata.obs["is_train"] = adata.obs["is_train"].astype(bool)
    adata.uns["method_code_version"] = check_r_version("Seurat")
    return adata
