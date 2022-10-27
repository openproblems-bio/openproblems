from ....data.tabula_muris_senis import load_tabula_muris_senis
from ....data.utils import filter_genes_cells
from ....tools.decorators import dataset
from .utils import generate_synthetic_dataset

import functools


def _tabula_muris_senis(alpha, test, n_obs):
    adata = load_tabula_muris_senis(
        test=test, organ_list=["lung"], method_list=["droplet"]
    )
    adata = adata[adata.obs["age"] == "30m"].copy()
    adata.obs["label"] = adata.obs["free_annotation"]

    merged_adata = generate_synthetic_dataset(
        adata, n_obs=n_obs, alpha=alpha, test=test
    )
    filter_genes_cells(merged_adata)
    return merged_adata


_tabula_muris_senis_dataset = functools.partial(
    dataset,
    data_url=load_tabula_muris_senis.metadata["data_url"],
    data_reference=load_tabula_muris_senis.metadata["data_reference"],
)


@_tabula_muris_senis_dataset(
    "Tabula muris senis (alpha=1)",
    dataset_summary="Mouse lung cells aggregated from single-cell"
    " (Dirichlet alpha=1)",
)
def tabula_muris_senis_alpha_1(test=False, n_obs=100):
    return _tabula_muris_senis(alpha=1, test=test, n_obs=n_obs)


@_tabula_muris_senis_dataset(
    "Tabula muris senis (alpha=5)",
    dataset_summary="Mouse lung cells aggregated from single-cell"
    " (Dirichlet alpha=5)",
)
def tabula_muris_senis_alpha_5(test=False, n_obs=100):
    return _tabula_muris_senis(alpha=5, test=test, n_obs=n_obs)


@_tabula_muris_senis_dataset(
    "Tabula muris senis (alpha=0.5)",
    dataset_summary="Mouse lung cells aggregated from single-cell"
    " (Dirichlet alpha=0.5)",
)
def tabula_muris_senis_alpha_0_5(test=False, n_obs=100):
    return _tabula_muris_senis(alpha=0.5, test=test, n_obs=n_obs)
