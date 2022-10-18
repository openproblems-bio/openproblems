from ...batch_integration_graph.methods.scalex import _scalex
from ...batch_integration_graph.methods.scalex import _scalex_method
from typing import Optional


@_scalex_method(method_name="SCALEX")
def scalex_full(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=False,
        compute_features=True,
        n_top_features=0,
    )
