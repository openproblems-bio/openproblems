from ....tools.decorators import method
from ....tools.normalize import preprocess_logCPM_1kHVG
from ....tools.utils import check_version


@method(
    method_name="Principle Component Analysis (PCA) (logCPM, 1kHVG)",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_url="https://www.tandfonline.com/doi/abs/10.1080/14786440109462720",
    paper_year=1901,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.decomposition.PCA.html",
    code_version=check_version("scikit-learn"),
)
def pca_logCPM_1kHVG(adata, test: bool = False):
    adata = preprocess_logCPM_1kHVG(adata)
    adata.obsm["X_emb"] = adata.obsm["X_input"][:, :2]
    return adata
