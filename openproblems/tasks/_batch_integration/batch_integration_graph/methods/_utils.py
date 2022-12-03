def hvg_batch(adata, batch_key, target_genes, adataOut):
    from scib.preprocessing import hvg_batch

    if adata.n_vars < 2000:
        return adata
    else:
        # uns and var get trampled
        uns = adata.uns.copy()
        var = adata.var.copy()
        adata = hvg_batch(
            adata,
            batch_key=batch_key,
            target_genes=target_genes,
            flavor="cell_ranger",
            adataOut=adataOut,
        )
        adata.var = var.loc[adata.var.index]
        adata.uns = uns
        return adata


def scale_batch(adata, batch_key):
    from scib.preprocessing import scale_batch

    # uns and var get trampled
    uns = adata.uns.copy()
    var = adata.var.copy()
    adata = scale_batch(adata, batch_key)
    adata.var = var.loc[adata.var_names]
    adata.uns = uns
    return adata
