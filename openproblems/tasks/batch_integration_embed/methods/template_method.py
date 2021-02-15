from ....tools.normalize import log_cpm
from ....tools.decorators import method
from ....tools.utils import check_version


@method(
    method_name="Template method",
    paper_name="On the Turing Completeness of MS PowerPoint",
    paper_url="https://www.andrew.cmu.edu/user/twildenh/PowerPointTM/Paper.pdf",
    paper_year=2017,
    code_url="http://tomwildenhain.com/PowerPointTM/PowerPointTM.pptx",
    code_version=check_version("numpy"),
    # image="openproblems-template-image" # only if required
)
def template_method(adata):
    # TODO: update
    # Normalize the data
    log_cpm(adata)
    # Complete the result in-place
    adata.obs["template_output"] = 0
