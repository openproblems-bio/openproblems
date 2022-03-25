from ....tools.decorators import method
from ....tools.utils import check_version
from .preprocessing import preprocess_scanpy


@method(
    method_name="Principle Component Analysis (PCA)",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_url="https://www.tandfonline.com/doi/abs/10.1080/14786440109462720",
    paper_year=1901,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.decomposition.PCA.html",
    code_version=check_version("scikit-learn"),
)
def pca(adata, test=False):
    preprocess_scanpy(adata)
    adata.obsm["X_emb"] = adata.obsm["X_input"][:, :2]
    return adata
