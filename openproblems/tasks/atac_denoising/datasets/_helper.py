import numpy as np


def _do_dropout(adata, seed, dropout_rate=0.3, cell_fraction=0.8):
    """
    Substract `dropout_rate` many atac reads from the positive
    peak counts.
    """
    np.random.seed(seed)

    n_cells = int(cell_fraction * adata.n_obs)
    affected_cells = np.random.choice(adata.n_obs, n_cells)
    X = adata.obsm["mode2"].copy()
    atac_subset = X[affected_cells, :]

    positive = atac_subset.data > 0
    values = np.asarray(atac_subset.data[positive].data, dtype=int)
    n_effects = int(dropout_rate * np.sum(values))
    choices = np.repeat(range(len(values)), values)
    affected_peaks = np.random.choice(choices, n_effects)
    idx, diff = np.unique(affected_peaks, return_counts=True)

    atac_subset.data[idx] -= diff
    X[affected_cells, :] = atac_subset
    adata.obsm["mode2_noisy"] = X
    return adata
