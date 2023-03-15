from ....tools.decorators import method
from ....tools.utils import check_version
from .xgboost import _xgboost
from typing import Optional

import functools

_scanvi_method = functools.partial(
    method,
    method_summary=(
        'scANVI or "single-cell ANnotation using Variational Inference" is a'
        " semi-supervised variant of the scVI(Lopez et al. 2018) algorithm. Like scVI,"
        " scANVI uses deep neural networks and stochastic optimization to model"
        " uncertainty caused by technical noise and bias in single - cell"
        " transcriptomics measurements. However, scANVI also leverages cell type labels"
        " in the generative modelling. In this approach, scANVI is used to predict the"
        " cell type labels of the unlabelled test data."
    ),
    paper_name=(
        "Probabilistic harmonization and annotation of single-cell transcriptomics data"
        " with deep generative models"
    ),
    paper_reference="xu2021probabilistic",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-pytorch",
)

_scanvi_scarches_method = functools.partial(
    method,
    method_summary=(
        'scArches+scANVI or "Single-cell architecture surgery" is a deep learning'
        " method for mapping new datasets onto a pre-existing reference model, using"
        " transfer learning and parameter optimization. It first uses scANVI to build a"
        " reference model from the training data, and then apply scArches to map the"
        " test data onto the reference model and make predictions."
    ),
    paper_name="Query to reference single-cell integration with transfer learning",
    paper_reference="lotfollahi2020query",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-pytorch",
)


def _hvg(adata, test=False, n_top_genes=2000):
    import scanpy as sc

    hvg_kwargs = dict(
        flavor="seurat_v3",
        inplace=False,
        n_top_genes=n_top_genes,
        batch_key="batch",
    )
    if test:
        hvg_kwargs["span"] = 0.8
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


def _scanvi(adata, test=False, n_hidden=None, n_latent=None, n_layers=None):
    import scvi

    if test:
        n_latent = n_latent or 10
        n_layers = n_layers or 1
        n_hidden = n_hidden or 32
    else:
        n_latent = n_latent or 30
        n_layers = n_layers or 2
        n_hidden = n_hidden or 128

    scanvi_labels = adata.obs["labels"].to_numpy()
    # test set labels masked
    scanvi_labels[~adata.obs["is_train"].to_numpy()] = "Unknown"
    adata.obs["scanvi_labels"] = scanvi_labels
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="scanvi_labels")
    scvi_model = scvi.model.SCVI(
        adata, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers
    )
    train_kwargs = dict(
        train_size=0.9,
        early_stopping=True,
    )
    if test:
        train_kwargs["max_epochs"] = 1
        train_kwargs["limit_train_batches"] = 10
        train_kwargs["limit_val_batches"] = 10
    scvi_model.train(**train_kwargs)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
    model.train(**train_kwargs)
    preds = model.predict(adata)
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return preds


def _scanvi_scarches(
    adata,
    test=False,
    n_hidden=None,
    n_latent=None,
    n_layers=None,
    prediction_method="scanvi",
):
    import numpy as np
    import scvi

    if test:
        n_latent = n_latent or 10
        n_layers = n_layers or 1
        n_hidden = n_hidden or 32
    else:
        n_latent = n_latent or 30
        n_layers = n_layers or 2
        n_hidden = n_hidden or 128

    unlabeled_category = "Unknown"

    # new obs labels to mask test set
    adata.obs["scanvi_labels"] = np.where(
        adata.obs["is_train"], adata.obs["labels"], unlabeled_category
    )
    adata_train = adata[adata.obs["is_train"]].copy()
    adata_test = adata[~adata.obs["is_train"]].copy()
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
    train_kwargs = dict(
        train_size=0.9,
        early_stopping=True,
    )
    if test:
        train_kwargs["max_epochs"] = 1
        train_kwargs["limit_train_batches"] = 10
        train_kwargs["limit_val_batches"] = 10
    scvi_model.train(**train_kwargs)
    model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category=unlabeled_category
    )
    model.train(**train_kwargs)

    query_model = scvi.model.SCANVI.load_query_data(adata_test, model)
    train_kwargs = dict(max_epochs=200, early_stopping=True)
    if test:
        train_kwargs["max_epochs"] = 1
        train_kwargs["limit_train_batches"] = 10
        train_kwargs["limit_val_batches"] = 10
    query_model.train(plan_kwargs=dict(weight_decay=0.0), **train_kwargs)

    if prediction_method == "scanvi":
        preds = _pred_scanvi(adata, query_model)
    elif prediction_method == "xgboost":
        preds = _pred_xgb(adata, query_model, test=test)

    return preds


def _pred_scanvi(adata, query_model):
    # this is temporary and won't be used
    adata.obs["scanvi_labels"] = "Unknown"
    preds = query_model.predict(adata)
    del adata.obs["scanvi_labels"]
    # predictions for train and test
    return preds


# note: could extend test option
def _pred_xgb(
    adata,
    query_model,
    test=False,
    num_round: Optional[int] = None,
):
    adata.obsm["X_emb"] = query_model.get_latent_representation(adata)
    adata = _xgboost(
        adata, test=test, obsm="X_emb", num_round=num_round, tree_method="hist"
    )
    return adata.obs["labels_pred"]


@_scanvi_method(method_name="scANVI (All genes)")
def scanvi_all_genes(adata, test=False):
    adata.obs["labels_pred"] = _scanvi(adata, test=test)
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_method(method_name="scANVI (Seurat v3 2000 HVG)")
def scanvi_hvg(adata, test=False):
    hvg_df = _hvg(adata, test)
    bdata = adata[:, hvg_df.highly_variable].copy()
    adata.obs["labels_pred"] = _scanvi(bdata, test=test)
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_scarches_method(method_name="scArches+scANVI (All genes)")
def scarches_scanvi_all_genes(adata, test=False):
    adata.obs["labels_pred"] = _scanvi_scarches(adata, test=test)
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_scarches_method(method_name="scArches+scANVI (Seurat v3 2000 HVG)")
def scarches_scanvi_hvg(adata, test=False):
    hvg_df = _hvg(adata, test)
    bdata = adata[:, hvg_df.highly_variable].copy()
    adata.obs["labels_pred"] = _scanvi_scarches(bdata, test=test)
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_scarches_method(method_name="scArches+scANVI+xgboost (All genes)")
def scarches_scanvi_xgb_all_genes(adata, test=False):
    adata.obs["labels_pred"] = _scanvi_scarches(
        adata, test=test, prediction_method="xgboost"
    )

    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_scarches_method(method_name="scArches+scANVI+xgboost (Seurat v3 2000 HVG)")
def scarches_scanvi_xgb_hvg(adata, test=False):
    hvg_df = _hvg(adata, test)
    bdata = adata[:, hvg_df.highly_variable].copy()
    adata.obs["labels_pred"] = _scanvi_scarches(
        bdata, test=test, prediction_method="xgboost"
    )

    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata
