import scanpy as sc
import scvi

def hvg_kwargs(n_top_genes=2000):
    return {
        "flavor": "seurat_v3",
        "inplace": False,
        "n_top_genes": n_top_genes,
        "batch_key": "batch",
    }

def hvg (adata, hvg_kwargs):
    try:
        return sc.pp.highly_variable_genes(adata[adata.obs["is_train"]], **hvg_kwargs)
    except ValueError:  # loess estimation can fail on small data with seurat_v3 flavor
        # in this case we try seurat flavor
        # and copy the data because it needs normalized counts
        # but later we need raw counts
        hvg_kwargs["flavor"] = "seurat"
        normdata = adata.copy()
        sc.pp.normalize_total(normdata, target_sum=1e4)
        sc.pp.log1p(normdata)
        return sc.pp.highly_variable_genes(
            normdata[normdata.obs["is_train"]], **hvg_kwargs
        )

def scanvi(adata, n_hidden=None, n_latent=None, n_layers=None, **train_kwargs):
    scanvi_labels = adata.obs["labels"].to_numpy()
    # test set labels masked
    scanvi_labels[~adata.obs["is_train"].to_numpy()] = "Unknown"
    adata.obs["scanvi_labels"] = scanvi_labels
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="scanvi_labels")
    scvi_model = scvi.model.SCVI(
        adata, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers
    )
    scvi_model.train(**train_kwargs)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train(**train_kwargs)
    preds = model.predict(adata)
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return preds
