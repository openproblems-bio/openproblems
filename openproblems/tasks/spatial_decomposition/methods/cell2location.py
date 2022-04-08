from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_version

# import cell2location


@method(
    method_name="Cell2Location",
    paper_name="Comprehensive mapping of tissue cell architecture via integrated single cell and spatial transcriptomics",  # noqa: E501
    paper_url="https://www.biorxiv.org/content/10.1101/2020.11.15.378125v1",
    paper_year=2020,
    code_url="https://github.com/BayraktarLab/cell2location",
    code_version=check_version("cell2location"),
)
def cell2location_log_cpm(adata, test=False):
    log_cpm(adata)
    # do something
    return
