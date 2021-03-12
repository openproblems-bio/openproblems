def _get_split(adata):
    uni = adata
    uni.X = uni.layers["logcounts"]
    uni.obsm["X_pca"] = uni.obsm["X_uni"]
    return (uni, adata)
