from ....tools.decorators import method
from ....tools.utils import check_version
from anndata import AnnData


@method(
    method_name="scvis",
    paper_name="Interpretable dimensionality reduction "
    "of single celltranscriptome data with deep generative models",
    paper_url="https://www.nature.com/articles/s41467-018-04368-5",
    paper_year=2018,
    code_url="https://bitbucket.org/jerry00/scvis-dev/",
    code_version=check_version("scvis"),
    image="openproblems-python-method-scvis",
)
def scvis(adata: AnnData) -> AnnData:
    return adata
