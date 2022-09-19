def _get_split(adata):
    uni = adata
    uni.obsm["X_pca"] = uni.obsm["X_uni_pca"]

    if "X_emb" not in adata.obsm:
        adata.obsm["X_emb"] = adata.obsm["X_pca"]

    return (uni, adata)
