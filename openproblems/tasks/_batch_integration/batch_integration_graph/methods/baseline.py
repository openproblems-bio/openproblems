from .....tools.decorators import method
from .....tools.utils import check_version

import numpy as np
import scipy.sparse


def _set_uns(adata):
    adata.uns["neighbors"] = adata.uns["uni"]
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    adata.uns["neighbors"]["distances_key"] = "distances"


@method(
    method_name="No Integration",
    paper_name="No Integration (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def no_integration(adata, test=False):
    adata.obsp["connectivities"] = adata.obsp["uni_connectivities"]
    adata.obsp["distances"] = adata.obsp["uni_distances"]
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


def _randomize_subgraph(distances, connectivities):
    idx = np.random.permutation(np.arange(distances.shape[0]))
    return distances[idx][:, idx], connectivities[idx][:, idx]


def _randomize_graph(distances, connectivities, partition=None):
    if partition is None:
        return _randomize_subgraph(distances, connectivities)
    else:
        distances_out = scipy.sparse.csr_matrix(distances.shape)
        connectivities_out = scipy.sparse.csr_matrix(connectivities.shape)
        for batch_name in np.unique(partition):
            idx = np.argwhere(partition == batch_name).flatten()
            distances_out[idx], connectivities_out[idx] = _randomize_subgraph(
                distances[idx], connectivities[idx]
            )
        return distances_out, connectivities_out


@method(
    method_name="Random Integration",
    paper_name="Random Integration (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_integration(adata, test=False):
    adata.obsp["distances"], adata.obsp["connectivities"] = _randomize_graph(
        adata.obsp["uni_distances"], adata.obsp["uni_connectivities"]
    )
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Integration by Celltype",
    paper_name="Random Integration by Celltype (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def celltype_random_integration(adata, test=False):
    adata.obsp["distances"], adata.obsp["connectivities"] = _randomize_graph(
        adata.obsp["uni_distances"],
        adata.obsp["uni_connectivities"],
        partition=adata.obs["labels"].to_numpy(),
    )
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Integration by Batch",
    paper_name="Random Integration by Batch (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def batch_random_integration(adata, test=False):
    adata.obsp["distances"], adata.obsp["connectivities"] = _randomize_graph(
        adata.obsp["uni_distances"],
        adata.obsp["uni_connectivities"],
        partition=adata.obs["batch"].to_numpy(),
    )
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
