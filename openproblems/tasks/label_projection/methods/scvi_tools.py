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
