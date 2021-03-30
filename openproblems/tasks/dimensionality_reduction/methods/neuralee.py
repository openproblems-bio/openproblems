from ....tools.decorators import method
from ....tools.utils import check_version
from anndata import AnnData


@method(
    method_name="",
    paper_name="",
    paper_url="",
    paper_year=2020,
    code_url="",
    code_version=check_version("neuralee"),
    image="openproblems-python-extras",
)
def neuralee(adata: AnnData) -> AnnData:
    return adata
