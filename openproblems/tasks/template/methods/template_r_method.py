from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

_template_method = r_function("template_r_method.R")


@method(
    method_name="Template method",
    paper_name="On the Turing Completeness of MS PowerPoint",
    paper_url="https://www.andrew.cmu.edu/user/twildenh/PowerPointTM/Paper.pdf",
    paper_year=2017,
    code_url="http://tomwildenhain.com/PowerPointTM/PowerPointTM.pptx",
    code_version=check_version("rpy2"),
    image="openproblems-r-base",
)
def template_method(adata):
    # TODO: update
    return _template_method(adata)
