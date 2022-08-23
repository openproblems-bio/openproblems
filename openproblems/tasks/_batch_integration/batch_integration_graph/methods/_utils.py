def hvg_batch(adata, batch_key, target_genes, adataOut):
    from scib.preprocessing import hvg_batch

    if adata.n_vars < 2000:
        return adata
    else:
        var = adata.var.copy()
        adata = hvg_batch(
            adata,
            batch_key=batch_key,
            target_genes=target_genes,
            flavor="cell_ranger",
            adataOut=adataOut,
        )
        adata.var = var.loc[adata.var.index]
        return adata


def scale_batch(adata, batch_key):
    from scib.preprocessing import scale_batch

    var = adata.var.copy()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata = scale_batch(adata, batch_key)
    adata.var = var.loc[adata.var_names]
    return adata
