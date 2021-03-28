from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="scGen",
    paper_name="scGen predicts single-cell perturbation responses",
    paper_url="https://www.nature.com/articles/s41592-019-0494-8",
    paper_year=2019,
    code_url="https://github.com/theislab/scgen",
    code_version=check_version("scgen"),
    image="openproblems-python37-scgen",  # only if required
)
def scgen_full_unscaled(adata):
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data

    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="scGen (hvg/unscaled)",
    paper_name="scGen predicts single-cell perturbation responses",
    paper_url="https://www.nature.com/articles/s41592-019-0494-8",
    paper_year=2019,
    code_url="https://github.com/theislab/scgen",
    code_version=check_version("scgen"),
    image="openproblems-python37-scgen",  # only if required
)
def scgen_hvg_unscaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="scGen (full/scaled)",
    paper_name="scGen predicts single-cell perturbation responses",
    paper_url="https://www.nature.com/articles/s41592-019-0494-8",
    paper_year=2019,
    code_url="https://github.com/theislab/scgen",
    code_version=check_version("scgen"),
    image="openproblems-python37-scgen",  # only if required
)
def scgen_full_scaled(adata):
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="scGen (hvg/scaled)",
    paper_name="scGen predicts single-cell perturbation responses",
    paper_url="https://www.nature.com/articles/s41592-019-0494-8",
    paper_year=2019,
    code_url="https://github.com/theislab/scgen",
    code_version=check_version("scgen"),
    image="openproblems-python37-scgen",  # only if required
)
def scgen_hvg_scaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata, umap=False)
    return adata
