from ...batch_integration_embed.metrics.utils import embedding_to_graph


def feature_to_embedding(adata):
    import scanpy as sc

    if adata.uns["is_baseline"] and "X_emb" in adata.obsm:
        # precomputed; do nothing
        return adata

    adata.obsm["X_emb"] = sc.pp.pca(adata.X)
    return adata


def feature_to_graph(adata):
    adata = feature_to_embedding(adata)
    adata = embedding_to_graph(adata)
    return adata
