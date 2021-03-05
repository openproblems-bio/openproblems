# from ....tools.normalize import log_cpm
from .....tools.decorators import method

# from ....tools.utils import check_version


@method(
    method_name="scanvi",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    # code_version=check_version("numpy"),
    image="openproblems-python-batch-integration" # only if required
)
def scanvi_emb(adata):
    from scIB.integration import runScanvi
    from scIB.preprocessing import reduce_data
    adata = runScanvi(adata, "batch", "labels")
    reduce_data(adata, use_emb="X_emb")
    return adata
