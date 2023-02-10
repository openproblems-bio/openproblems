def feature_to_embedding(adata):
    import scanpy as sc

    if adata.uns["is_baseline"] and "X_emb" in adata.obsm:
        # precomputed; do nothing
        return adata

    adata.obsm["X_emb"] = sc.pp.pca(adata.X)
    return adata


def embedding_to_graph(adata):
    import scanpy as sc

    if adata.uns["is_baseline"] and "neighbors" in adata.uns:
        # precomputed; do nothing
        return adata

    sc.pp.neighbors(adata, use_rep="X_emb")
    return adata


def feature_to_graph(adata):
    adata = feature_to_embedding(adata)
    adata = embedding_to_graph(adata)
    return adata
