from ....tools.decorators import method
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version
from .preprocessing import preprocess_scanpy


def _phate(adata):
    from phate import PHATE

    phate_op = PHATE(verbose=False, n_jobs=-1)
    adata.obsm["X_emb"] = phate_op.fit_transform(adata.X)
    return adata


@method(
    method_name="PHATE (default pre-processing)",
    paper_name="Visualizing Transitions and Structure for Biological Data Exploration",
    paper_url="https://www.nature.com/articles/s41587-019-0336-3",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    code_version=check_version("phate"),
    image="openproblems-python-extras",
)
def phate_default(adata):
    sqrt_cpm(adata)
    return _phate(adata)


@method(
    method_name="PHATE (scanpy pre-processing)",
    paper_name="Visualizing Transitions and Structure for Biological Data Exploration",
    paper_url="https://www.nature.com/articles/s41587-019-0336-3",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    code_version=check_version("phate"),
    image="openproblems-python-extras",
)
def phate_scanpy(adata):
    preprocess_scanpy(adata)
    return _phate(adata)
