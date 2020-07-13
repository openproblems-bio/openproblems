from sklearn.decomposition import TruncatedSVD
import scipy.spatial


def procrustes(adata):
    X_pca = TruncatedSVD(100).fit_transform(adata.X)
    Y_pca = TruncatedSVD(100).fit_transform(adata.obsm["mode2"])
    X_proc, Y_proc, _ = scipy.spatial.procrustes(X_pca, Y_pca)
    adata.obsm["aligned"] = X_proc
    adata.obsm["mode2_aligned"] = Y_proc
