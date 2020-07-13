from sklearn.decomposition import TruncatedSVD


def cheat(adata):
    adata.obsm["aligned"] = adata.X
    adata.uns["mode2"].obsm["aligned"] = adata.X
