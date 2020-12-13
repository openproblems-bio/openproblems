from ....tools.decorators import method
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version

from phate import PHATE


@method(
    method_name="PHATE",
    paper_name="Visualizing Transitions and Structure for Biological Data Exploration",
    paper_url="https://www.nature.com/articles/s41587-019-0336-3",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    code_version=check_version("phate"),
)
def phate(adata):
    sqrt_cpm(adata)
    phate_op = PHATE(verbose=False, n_jobs=-1)
    adata.obsm["X_emb"] = phate_op.fit_transform(adata.X)
    return adata
