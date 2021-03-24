# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Combat",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_full_unscaled(adata):
    from scIB.integration import runCombat
    from scIB.preprocessing import reduce_data

    adata = runCombat(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata


@method(
    method_name="Combat (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("combatpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_hvg_unscaled(adata):
    from ._hvg import hvg_batch
    from scIB.integration import runCombat
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runCombat(adata, "batch")
    reduce_data(adata)
    return adata


@method(
    method_name="Combat (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("combatpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_hvg_scaled(adata):
    from ._hvg import hvg_batch
    from ._hvg import scale_batch
    from scIB.integration import runCombat
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runCombat(adata, "batch")
    reduce_data(adata)
    return adata


@method(
    method_name="Combat (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("combatpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_full_scaled(adata):
    from ._hvg import scale_batch
    from scIB.integration import runCombat
    from scIB.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = runCombat(adata, "batch")
    reduce_data(adata)
    return adata
