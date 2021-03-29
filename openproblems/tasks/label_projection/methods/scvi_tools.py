from ....tools.decorators import method
from ....tools.utils import check_version


def _scanvi(adata):
    import scvi

    scanvi_labels = adata.obs["labels"].to_numpy()
    # test set labels masked
    scanvi_labels[~adata.obs["is_train"].to_numpy()] = "Unknown"
    adata.obs["scanvi_labels"] = scanvi_labels
    scvi.data.setup_anndata(adata, batch_key="batch", labels_key="scanvi_labels")
    scvi_model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
    scvi_model.train(train_size=1.0)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train(train_size=1.0)
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return model.predict(adata)


def _scanvi_scarches(adata):
    import scvi

    # new obs labels to mask test set
    adata_train = adata[adata.obs["is_train"]].copy()
    adata_train.obs["scanvi_labels"] = adata_train.obs["labels"].copy()
    adata_test = adata[~adata.obs["is_train"]].copy()
    adata_test.obs["scanvi_labels"] = "Unknown"
    scvi.data.setup_anndata(adata_train, batch_key="batch", labels_key="scanvi_labels")

    # specific scArches parameters
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        n_latent=30,
    )
    scvi_model = scvi.model.SCVI(adata, **arches_params)
    scvi_model.train(train_size=1.0)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train(train_size=1.0)

    query_model = scvi.model.SCANVI.load_query_data(adata_test, model)
    query_model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))

    # predictions for train and test
    return query_model.predict(adata)


@method(
    method_name="scANVI (All genes)",
    paper_name="Probabilistic harmonization and annotation of single-cell"
    " transcriptomics data with deep generative models.",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-scvi",
)
def scanvi_all_genes(adata):
    adata.obs["labels_pred"] = _scanvi(adata)
    return adata


@method(
    method_name="scANVI (Seurat v3 2000 HVG)",
    paper_name="Probabilistic harmonization and annotation of single-cell"
    " transcriptomics data with deep generative models.",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-scvi",
)
def scanvi_hvg(adata):
    import scanpy as sc

    hvg_df = sc.pp.highly_variable_genes(
        adata[adata.obs["is_train"]],
        flavor="seurat_v3",
        inplace=False,
        n_top_genes=2000,
        batch_key="batch",
    )
    bdata = adata[:, hvg_df.highly_variable].copy()
    adata.obs["labels_pred"] = _scanvi(bdata)
    return adata


@method(
    method_name="scArches+scANVI (All genes)",
    paper_name="Query to reference single-cell integration with transfer learning.",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.07.16.205997v1",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-scvi",
)
def scarches_scanvi_all_genes(adata):
    adata.obs["labels_pred"] = _scanvi_scarches(adata)
    return adata


@method(
    method_name="scArches+scANVI (Seurat v3 2000 HVG)",
    paper_name="Query to reference single-cell integration with transfer learning.",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.07.16.205997v1",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-scvi",
)
def scarches_scanvi_hvg(adata):
    import scanpy as sc

    hvg_df = sc.pp.highly_variable_genes(
        adata[adata.obs["is_train"]],
        flavor="seurat_v3",
        inplace=False,
        n_top_genes=2000,
        batch_key="batch",
    )
    bdata = adata[:, hvg_df.highly_variable].copy()
    adata.obs["labels_pred"] = _scanvi_scarches(bdata)
    return adata
