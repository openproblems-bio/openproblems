from ....data.pancreas import get_pancreas_integer
from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset
from .utils import generate_synthetic_dataset

import scanpy as sc


@dataset(
    "Pancreas (alpha=1)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_alpha_1(test=False, n_obs=100):
    adata = load_pancreas(test=test)
    adata = get_pancreas_integer(adata)
    sc.pp.filter_genes(adata, min_counts=10)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=1)
    return merged_adata


@dataset(
    "Pancreas (alpha=5)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_alpha_5(test=False, n_obs=100):
    adata = load_pancreas(test=test)
    adata = get_pancreas_integer(adata)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=5)
    return merged_adata


@dataset(
    "Pancreas (alpha=0.5)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_alpha_0_5(test=False, n_obs=100):
    adata = load_pancreas(test=test)
    adata = get_pancreas_integer(adata)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(adata, n_obs=n_obs, alpha=0.5)
    return merged_adata
