def hvg_batch(adata, batch_key, target_genes, adataOut):
    from scib.preprocessing import hvg_batch

    if adata.n_vars < 2000:
        return adata
    else:
        return hvg_batch(
            adata,
            batch_key=batch_key,
            target_genes=target_genes,
            flavor="cell_ranger",
            adataOut=adataOut,
        )


def scale_batch(adata, batch_key):
    from scib.preprocessing import scale_batch

    adata.strings_to_categoricals()
    tmp = scale_batch(adata, batch_key)
    tmp.var = tmp.var.iloc[:, :-6]
    return tmp
