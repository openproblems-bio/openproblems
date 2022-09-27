from .....tools.decorators import method
from .....tools.utils import check_version

import functools

_combat_method = functools.partial(
    method,
    paper_name="Adjusting batch effects in microarray expression data using "
    "empirical Bayes methods",
    paper_url="https://academic.oup.com/biostatistics/article/8/1/118/252073",
    paper_year=2007,
    code_url="https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html",
    image="openproblems-python-batch-integration",
)


def _combat(adata):
    from openproblems.tools.normalize import log_scran_pooling
    from scib.integration import runCombat
    from scib.preprocessing import reduce_data

    adata = runCombat(adata, "batch")
    reduce_data(adata, umap=False)
    adata.uns["method_code_version"] = check_version("scanpy")
    return adata


@_combat_method(method_name="Combat (full/unscaled)")
def combat_full_unscaled(adata, test=False):
    return _combat(adata)


@_combat_method(method_name="Combat (hvg/unscaled)")
def combat_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _combat(adata)


@_combat_method(method_name="Combat (hvg/scaled)")
def combat_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _combat(adata)


@_combat_method(method_name="Combat (full/scaled)")
def combat_full_scaled(adata, test=False):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _combat(adata)
