from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_version

# from scvi.external import SpatialStereoscope
# from scvi.external import RNAStereoscope


@method(
    method_name="Stereoscope",
    paper_name="Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography",  # noqa: E501
    paper_url="https://www.nature.com/articles/s42003-020-01247-y",
    paper_year=2020,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi-tools"),
)
def stereoscope_log_cpm(adata):
    log_cpm(adata)
    # do something
    return
