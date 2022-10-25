from ....data.pancreas import load_pancreas
from ....data.utils import filter_genes_cells
from ....tools.decorators import dataset
from .utils import generate_synthetic_dataset

import functools
import scanpy as sc

_pancreas_dataset = functools.partial(
    dataset,
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
)
_DATASET_SUMMARY = (
    "Human pancreas cells aggregated from single-cell (Dirichlet alpha={})"
)


def _pancreas_synthetic(alpha: float, test: bool = False, n_obs: int = 100):
    adata = load_pancreas(test=test, keep_techs=["inDrop3"])
    sc.pp.filter_genes(adata, min_counts=10)
    adata.obs["label"] = adata.obs["celltype"]

    merged_adata = generate_synthetic_dataset(
        adata, n_obs=n_obs, alpha=alpha, test=test
    )
    filter_genes_cells(merged_adata)
    return merged_adata


@_pancreas_dataset(
    "Pancreas (alpha=1)",
    dataset_summary=_DATASET_SUMMARY.format(1),
)
def pancreas_alpha_1(test=False, n_obs=100):
    return _pancreas_synthetic(test=test, n_obs=n_obs, alpha=1)


@_pancreas_dataset(
    "Pancreas (alpha=5)",
    dataset_summary=_DATASET_SUMMARY.format(5),
)
def pancreas_alpha_5(test=False, n_obs=100):
    return _pancreas_synthetic(test=test, n_obs=n_obs, alpha=5)


@_pancreas_dataset(
    "Pancreas (alpha=0.5)",
    dataset_summary=_DATASET_SUMMARY.format(0.5),
)
def pancreas_alpha_0_5(test=False, n_obs=100):
    return _pancreas_synthetic(test=test, n_obs=n_obs, alpha=0.5)
