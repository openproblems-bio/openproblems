import numpy as np


def add_celltype_noise(adata, noise_prob):
    """Inject random celltype noise in the dataset .
    This is done by permuting a fraction of the celltypes in the training set.
    By adding different levels of celltype noise metrics can be evaluated to show
    generalization trends from training data even if ground truth is uncertain.
    Parameters
    -------
    adata : AnnData
        A dataset with the required fields for the celltype_projection task.
    noise_prob : Float
        The probability of celltype noise in the training data.
    Returns
    -------
    new_adata : AnnData
        Dataset where training celltypes have been permuted by specified probability.
    """

    old_celltype = adata.obs["celltype"].pipe(np.array)
    old_celltype_train = old_celltype[adata.obs["is_train"]].copy()
    new_celltype_train = old_celltype_train.copy()

    celltype_names = np.unique(new_celltype_train)

    n_celltype = celltype_names.shape[0]

    reassign_probs = (noise_prob / (n_celltype - 1)) * np.ones((n_celltype, n_celltype))

    np.fill_diagonal(reassign_probs, 1 - noise_prob)

    for k, celltype in enumerate(celltype_names):
        celltype_indices = np.where(old_celltype_train == celltype)[0]
        new_celltype_train[celltype_indices] = np.random.choice(
            celltype_names, celltype_indices.shape[0], p=reassign_probs[:, k]
        )

    adata.obs.loc[adata.obs["is_train"], "celltype"] = new_celltype_train

    return adata
