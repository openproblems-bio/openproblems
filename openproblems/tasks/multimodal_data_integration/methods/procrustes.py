from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_version


@method(
    method_name="Procrustes",
    paper_name="Generalized Procrustes analysis",
    paper_url="https://link.springer.com/content/pdf/10.1007/BF02291478.pdf",
    paper_year=1975,
    code_url="https://docs.scipy.org/doc/scipy/reference/generated/"
    "scipy.spatial.procrustes.html",
)
def procrustes(adata, test=False, n_svd=None):
    import scipy.spatial
    import sklearn.decomposition

    if test:
        n_svd = n_svd or 20
    else:  # pragma: no cover
        n_svd = n_svd or 100
    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    adata = log_cpm(adata)
    adata = log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])
    X_proc, Y_proc, _ = scipy.spatial.procrustes(X_pca, Y_pca)
    adata.obsm["aligned"] = X_proc
    adata.obsm["mode2_aligned"] = Y_proc

    adata.uns["method_code_version"] = check_version("scipy")
    return adata
