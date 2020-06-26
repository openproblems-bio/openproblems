from sklearn.decomposition import TruncatedSVD


def cheat(adata, n_svd=150):
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    adata.obsm["aligned"] = X_pca
    adata.uns["mode2"].obsm["aligned"] = X_pca
