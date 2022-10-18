from ...batch_integration_graph.methods.scalex import _scalex
from ...batch_integration_graph.methods.scalex import _scalex_method
from typing import Optional


@_scalex_method(method_name="SCALEX (full)")
def scalex_full(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=False,
        compute_features=False,
        n_top_features=0,
    )


@_scalex_method(method_name="SCALEX (hvg)")
def scalex_hvg(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=False,
        compute_features=False,
        n_top_features=2000,
    )
