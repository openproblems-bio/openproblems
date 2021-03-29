from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_version

# import rctd, this is an R tool so need to have also R file


@method(
    method_name="Rctd",
    paper_name="Robust decomposition of cell type mixtures in spatial transcriptomics",
    paper_url="https://www.nature.com/articles/s41587-021-00830-w",
    paper_year=2020,
    code_url="https://github.com/dmcable/RCTD",
    code_version=check_version("rctd"),
)
def stereoscope_log_cpm(adata):
    log_cpm(adata)
    # do something
    return