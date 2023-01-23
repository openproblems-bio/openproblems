def _get_split(adata):
    uni = adata
    uni.obsm["X_pca"] = uni.obsm["X_uni_pca"]
    uni.X = uni.layers["log_normalized"]
    return (uni, adata)
