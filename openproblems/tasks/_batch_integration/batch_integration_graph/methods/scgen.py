# from ....tools.normalize import log_cpm
from .....tools.decorators import method

# from ....tools.utils import check_version

@method(
    method_name="scGen",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    # code_version=check_version("numpy"),
    image="openproblems-python-batch-integration" # only if required
)
def scgen_full(adata):
    from scIB.integration import runScGen
    from scIB.preprocessing import reduce_data

    adata = runScGen(adata, "batch", "labels")
    reduce_data(adata)
    return adata
