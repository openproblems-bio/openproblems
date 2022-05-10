# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

@method(
    method_name="BBKNN",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def bbknn_full_unscaled(adata, test=False):
    from scib.integration import combat
    from scib.preprocessing import reduce_data

    adata = combat(adata, "batch")
    reduce_data(adata, umap=False)
    # Complete the result in-place
    return adata
