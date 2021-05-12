from ....tools.decorators import method
from ....tools.utils import check_version
from .preprocessing import preprocess_scanpy


@method(
    method_name="densMAP",
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    code_version=check_version("umap-learn"),
    image="openproblems-python-extras",
)
def densmap(adata):
    from umap import UMAP

    preprocess_scanpy(adata)
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(adata.X)
    return adata
