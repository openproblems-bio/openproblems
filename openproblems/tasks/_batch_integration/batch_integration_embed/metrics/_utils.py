def _get_split(adata):
    uni = adata
    uni.X = uni.layers["logcounts"]
    uni.obsm["X_pca"] = uni.obsm["X_uni"]

    if "X_emb" not in adata.obsm:
        adata.obsm["X_emb"] = adata.obsm["X_pca"]

    return (uni, adata)
