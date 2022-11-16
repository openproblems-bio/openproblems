import scib


def hvg_batch(adata, batch_key, n_hvg):
    """
    Compute highly variable genes by batch
    """
    if n_hvg > adata.n_vars or n_hvg == 0:
        return adata.var_names.tolist()
    return scib.pp.hvg_batch(
        adata,
        batch_key=batch_key,
        target_genes=n_hvg,
        adataOut=False
    )
