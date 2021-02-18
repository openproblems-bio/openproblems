# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from ....tools.utils import check_version
from scIB.integration import runBBKNN
from scIB.preprocessing import hvg_batch, scale_batch


@method(
    method_name="BBKNN",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="",
    code_version=check_version("bbknn"),
    # image="openproblems-template-image" # only if required
)
def bbknn_full_unscaled(adata):
    # Normalize the data

    runBBKNN(adata, "batch")
    # Complete the result in-place
    adata.obs["template_output"] = 0
    return adata
