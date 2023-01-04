from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version


@method(
    method_name="Principle Component Analysis (PCA) (logCPM, 1kHVG)",
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_reference="pearson1901pca",
    paper_year=1901,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.decomposition.PCA.html",
)
def pca_logCPM_1kHVG(adata, test: bool = False):
    import scanpy as sc

    adata = log_cpm_hvg(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    adata.obsm["X_emb"] = adata.obsm["X_pca"][:, :2]
    adata.uns["method_code_version"] = check_version("scikit-learn")
    return adata
