from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import scipy.sparse


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
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Integration",
    paper_name="Random Integration (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_integration(adata, test=False):
    distances, connectivities = (
        adata.obsp["uni_distances"],
        adata.obsp["uni_connectivities"],
    )
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
    adata.obsp["distances"] = scipy.sparse.coo_matrix(
        (distances_data, (row, col)), shape=distances.shape
    )
    connectivities_data = np.random.permutation(connectivities.data)[: len(row)]
    adata.obsp["connectivities"] = scipy.sparse.coo_matrix(
        (connectivities_data, (row, col)), shape=connectivities.shape
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
