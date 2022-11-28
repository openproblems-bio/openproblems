def _get_split(adata):
    uni = adata
    uni.obsm["X_pca"] = uni.obsm["X_uni_pca"]
    return (uni, adata)
