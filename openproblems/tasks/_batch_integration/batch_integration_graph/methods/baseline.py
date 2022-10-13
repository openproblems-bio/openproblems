from .....tools.decorators import method
from .....tools.utils import check_version

import numpy as np
import scipy.sparse


def _set_uns(adata):
    adata.uns["neighbors"] = adata.uns["uni_neighbors"]
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
    row = np.random.choice(
        np.arange(distances.shape[0]), len(distances.data), replace=True
    )
    col = np.random.choice(
        np.arange(distances.shape[0]), len(distances.data), replace=True
    )
    while np.any(np.any(row == col)):
        # eliminate self-loops
        row[row == col] = np.random.choice(
            np.arange(distances.shape[0]), np.sum(row == col)
        )
    # eliminate duplicates
    row, col = zip(*set(zip(row, col)))
    distances_data = np.random.permutation(distances.data)[: len(row)]
    distances_out = scipy.sparse.coo_matrix(
        (distances_data, (row, col)), shape=distances.shape
    ).tocsr()
    connectivities_data = np.random.permutation(connectivities.data)[: len(row)]
    connectivities_out = scipy.sparse.coo_matrix(
        (connectivities_data, (row, col)), shape=connectivities.shape
    ).tocsr()
    return distances_out, connectivities_out


def _randomize_graph(distances, connectivities, batch=None):
    if batch is None:
        return _randomize_subgraph(distances, connectivities)
    else:
        distances_out = scipy.sparse.csr_matrix(distances.shape)
        connectivities_out = scipy.sparse.csr_matrix(connectivities.shape)
        for batch_name in np.unique(batch):
            idx = np.argwhere(batch == batch_name).flatten()
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
def celltype_integration(adata, test=False):
    adata.obsp["distances"], adata.obsp["connectivities"] = _randomize_graph(
        adata.obsp["uni_distances"],
        adata.obsp["uni_connectivities"],
        batch=adata.obs["batch"].to_numpy(),
    )
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
