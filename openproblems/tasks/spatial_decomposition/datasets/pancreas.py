from ....data.pancreas import load_pancreas
from ....data.utils import filter_genes_cells
from ....tools.decorators import dataset
from .utils import generate_synthetic_dataset

import scanpy as sc


@dataset(
    "Pancreas (alpha=1)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreas cells aggregated from single-cell"
    " (Dirichlet alpha=1)",
)
def pancreas_alpha_1(test=False, n_obs=100):
    adata = load_pancreas(test=test, integer_only=True)
    sc.pp.filter_genes(adata, min_counts=10)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=1)
    filter_genes_cells(merged_adata)
    return merged_adata


@dataset(
    "Pancreas (alpha=5)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreas cells aggregated from single-cell"
    " (Dirichlet alpha=5)",
)
def pancreas_alpha_5(test=False, n_obs=100):
    adata = load_pancreas(test=test, integer_only=True)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=5)
    filter_genes_cells(merged_adata)
    return merged_adata


@dataset(
    "Pancreas (alpha=0.5)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreas cells aggregated from single-cell"
    " (Dirichlet alpha=0.5)",
)
def pancreas_alpha_0_5(test=False, n_obs=100):
    adata = load_pancreas(test=test, integer_only=True)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=0.5)
    filter_genes_cells(merged_adata)
    return merged_adata
