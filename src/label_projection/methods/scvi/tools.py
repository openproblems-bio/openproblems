import scanpy as sc
import scvi


def hvg (adata, **hvg_kwargs):
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
    scanvi_labels = adata.obs["celltype"].to_numpy()
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
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return model.predict(adata)

def scanvi_scarches(adata, n_hidden=None, n_latent=None, n_layers=None, train_kwargs={}):
    model_train_kwargs = train_kwargs['model_train_kwargs']
    query_model_train_kwargs = train_kwargs['query_model_train_kwargs']
    # new obs labels to mask test set
    adata_train = adata[adata.obs["is_train"]].copy()
    adata_train.obs["scanvi_labels"] = adata_train.obs["celltype"].copy()
    adata_test = adata[~adata.obs["is_train"]].copy()
    adata_test.obs["scanvi_labels"] = "Unknown"
    scvi.model.SCVI.setup_anndata(
        adata_train, batch_key="batch", labels_key="scanvi_labels"
    )

    # specific scArches parameters
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_hidden=n_hidden,
        n_layers=n_layers,
        n_latent=n_latent,
    )
    scvi_model = scvi.model.SCVI(adata_train, **arches_params)

    scvi_model.train(**model_train_kwargs)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train(**model_train_kwargs)

    query_model = scvi.model.SCANVI.load_query_data(adata_test, model)
    query_model.train(plan_kwargs=dict(weight_decay=0.0), **query_model_train_kwargs)

    # this is temporary and won't be used
    adata.obs["scanvi_labels"] = "Unknown"
    preds = query_model.predict(adata)
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return preds
