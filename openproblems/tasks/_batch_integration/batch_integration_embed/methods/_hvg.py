def hvg_batch(adata, batch_key, target_genes, adataOut):
    from scIB.preprocessing import hvg_batch

    if adata.n_vars < 2000:
        return adata
    else:
        hvg_batch(adata, batch_key, target_genes, adataOut)
