from sklearn.decomposition import TruncatedSVD


def cheat(adata):
    adata.obsm["aligned"] = adata.X
    adata.obsm["mode2_aligned"] = adata.X
