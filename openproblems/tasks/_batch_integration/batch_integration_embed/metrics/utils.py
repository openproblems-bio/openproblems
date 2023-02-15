def embedding_to_graph(adata):
    import scanpy as sc

    if adata.uns["is_baseline"] and "neighbors" in adata.uns:
        # precomputed; do nothing
        return adata

    sc.pp.neighbors(adata, use_rep="X_emb")
    return adata


def get_split(adata):
    uni = adata
    uni.obsm["X_pca"] = uni.obsm["X_uni_pca"]
    uni.X = uni.layers["log_normalized"]
    return (uni, adata)
