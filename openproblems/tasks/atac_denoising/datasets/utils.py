import anndata
import numpy as np


def do_dropout(adata, seed, dropout_rate=0.2, cell_fraction=0.8):
    """Set counts to 0 for `dropout_rate` many peaks."""
    np.random.seed(seed)

    n_cells = int(cell_fraction * adata.n_obs)
    affected_cells = np.random.choice(adata.n_obs, n_cells)

    X = adata.obsm["mode2"].copy()
    atac_subset = X[affected_cells, :]
    n_effects = int(dropout_rate * len(atac_subset.data))
    dropouts = np.random.choice(len(atac_subset.data), n_effects)

    atac_subset.data[dropouts] = 0
    X[affected_cells, :].data = atac_subset
    X.eliminate_zeros()
    adata.obsm["mode2_noisy"] = X

    return adata


def split_data(
    adata: anndata.AnnData, train_frac: float = 0.9, seed: int = 0
) -> anndata.AnnData:
    """Split data using molecular cross-validation."""
    import molecular_cross_validation.util
    import scipy

    random_state = np.random.RandomState(seed)

    X = adata.obsm["mode2"]

    if scipy.sparse.issparse(X):
        X = np.array(X.todense())
    if np.allclose(X, X.astype(np.int)):
        X = X.astype(np.int)
    else:
        raise TypeError("Molecular cross-validation requires integer count data.")

    X_train, X_test = molecular_cross_validation.util.split_molecules(
        X, 0.9, 0.0, random_state
    )
    adata.obsm["mode2_noisy"] = X_train
    return adata
