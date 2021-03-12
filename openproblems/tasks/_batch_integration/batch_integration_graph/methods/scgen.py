# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="scGen",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    code_version=check_version("scgen"),
    image="openproblems-python-scgen",  # only if required
)
def scgen_full_unscaled(adata):
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data

    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata)
    return adata


@method(
    method_name="scGen (hvg/unscaled)",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    code_version=check_version("scgen"),
    image="openproblems-python-scgen",  # only if required
)
def scgen_hvg_unscaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata)
    return adata


@method(
    method_name="scGen (full/scaled)",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    code_version=check_version("scgen"),
    image="openproblems-python-scgen",  # only if required
)
def scgen_full_scaled(adata):
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata)
    return adata


@method(
    method_name="scGen (hvg/scaled)",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    code_version=check_version("scgen"),
    image="openproblems-python-scgen",  # only if required
)
def scgen_hvg_scaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata)
    return adata
