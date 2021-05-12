from ....tools.decorators import method
from ....tools.utils import check_version
from .preprocessing import preprocess_scanpy
from anndata import AnnData


@method(
    method_name="Ivis (CPU)",
    paper_name="Structure-preserving visualisation of "
    "high dimensional single-cell datasets",
    paper_url="https://www.nature.com/articles/s41598-019-45301-0",
    paper_year=2019,
    code_url="https://github.com/beringresearch/ivis",
    code_version=check_version("ivis"),
    image="openproblems-python-extras",
)
def ivis(adata: AnnData) -> AnnData:
    from ivis import Ivis

    preprocess_scanpy(adata)

    # parameters taken from:
    # https://bering-ivis.readthedocs.io/en/latest/
    # scanpy_singlecell.html#reducing-dimensionality-using-ivis
    ivis = Ivis(
        k=15,
        model="maaten",
        n_epochs_without_progress=5,
        verbose=0,
        embedding_dims=2,
    )
    adata.obsm["X_emb"] = ivis.fit_transform(adata.obsm["X_pca"])

    return adata
