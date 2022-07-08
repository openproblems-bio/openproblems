from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import scanpy as sc


@method(
    method_name="Principle Component Analysis (PCA) (logCPM, 1kHVG)",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_url="https://www.tandfonline.com/doi/abs/10.1080/14786440109462720",
    paper_year=1901,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.decomposition.PCA.html",
)
def pca_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    adata.obsm["X_emb"] = adata.obsm["X_pca"][:, :2]
    adata.uns["method_code_version"] = check_version("scikit-learn")
    return adata
