from ....tools.decorators import method
from ....tools.utils import check_version

import scanpy as sc
import scvi


def _scanvi(adata):
    adata_train = adata[adata.obs["is_train"]].copy()
    scvi.data.setup_anndata(adata_train, batch_key="batch", labels_key="labels")
    scvi_model = scvi.model.SCVI(adata_train, n_latent=30, n_layers=2, train_size=1.0)
    scvi_model.train()
    # all cells treated as labeled
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train()
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
    bdata = adata.copy()
    bdata = sc.pp.highly_variable_genes(
        bdata, flavor="seurat_v3", subset=True, n_top_genes=2000
    )
    adata.obs["labels_pred"] = _scanvi(bdata)
    return adata
