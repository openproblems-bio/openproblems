from sklearn.decomposition import TruncatedSVD


def cheat(adata, n_svd=300):
    if min(adata.X.shape) <= n_svd:
        n_svd = min(adata.X.shape) - 1
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    adata.obsm["aligned"] = X_pca
    adata.obsm["mode2_aligned"] = X_pca
