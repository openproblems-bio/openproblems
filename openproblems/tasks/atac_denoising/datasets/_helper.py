import numpy as np


def _do_dropout(adata, seed, dropout_rate=0.2, cell_fraction=0.8):
    """Set counts to 0 for `dropout_rate` many peaks.
    """
    np.random.seed(seed)

    n_cells = int(cell_fraction * adata.n_obs)
    affected_cells = np.random.choice(adata.n_obs, n_cells)

    X = adata.obsm["mode2"].copy()
    atac_subset = X[affected_cells, :].data
    n_effects = int(dropout_rate * len(atac_subset))
    dropouts = np.random.choice(len(atac_subset), n_effects)

    atac_subset[dropouts] = 0
    X[affected_cells, :].data = atac_subset
    X.eliminate_zeros()
    adata.obsm["mode2_noisy"] = X

    return adata
