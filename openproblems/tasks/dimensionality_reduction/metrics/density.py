from ....tools.decorators import metric
from anndata import AnnData


@metric("density preservation", maximize=True)
def density_preservation(_adata: AnnData) -> float:
    return 0.0
