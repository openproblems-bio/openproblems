# from ....tools.normalize import log_cpm
from .....tools.decorators import method

# from ....tools.utils import check_version
from scIB.integration import runScanorama
from scIB.preprocessing import reduce_data


@method(
    method_name="Scanorama",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    # code_version=check_version("numpy"),
    # image="openproblems-template-image" # only if required
)
def scanorama_full(adata):
    runScanorama(adata, "batch")
    reduce_data(adata)
    return adata


def scanorama_emb(adata):
    runScanorama(adata, "batch")
    reduce_data(adata, use_rep="X_emb")
    return adata
