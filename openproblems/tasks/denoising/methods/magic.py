import magic

from ....tools.decorators import method
from ....tools.utils import check_version


@method(
    method_name="MAGIC",
    paper_name="Recovering Gene Interactions from Single-Cell Data "
    "Using Data Diffusion",
    paper_url="https://www.cell.com/cell/abstract/S0092-8674(18)30724-4",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    code_version=check_version("magic-impute"),
    image="openproblems-python-extras",
)
def no_denoising(adata):
    adata.obsm["denoised"] = magic.MAGIC().fit_transform(
        adata.obsm["train"], genes="all_genes"
    )
    return adata
