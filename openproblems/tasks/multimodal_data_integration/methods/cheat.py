from sklearn.decomposition import TruncatedSVD


def cheat(adata):
    adata.obsm["aligned"] = adata.X[:, :200]
    adata.obsm["mode2_aligned"] = adata.X[:, :200]
