import numpy as np


def add_label_noise(adata, noise_prob):
    """Inject random label noise in the dataset .
    This is done by permuting a fraction of the labels in the training set.
    By adding different levels of label noise metrics can be evaluated to show
    generalization trends from training data even if ground truth is uncertain.
    Parameters
    -------
    adata : AnnData
        A dataset with the required fields for the label_projection task.
    noise_prob : Float
        The probability of label noise in the training data.
    Returns
    -------
    new_adata : AnnData
        Dataset where training labels have been permuted by specified probability.
    """

    old_labels = adata.obs["labels"].pipe(np.array)
    old_labels_train = old_labels[adata.obs["is_train"]].copy()
    new_labels_train = old_labels_train.copy()

    label_names = np.unique(new_labels_train)

    n_labels = label_names.shape[0]

    reassign_probs = (noise_prob / (n_labels - 1)) * np.ones((n_labels, n_labels))

    np.fill_diagonal(reassign_probs, 1 - noise_prob)

    for k, label in enumerate(label_names):
        label_indices = np.where(old_labels_train == label)[0]
        new_labels_train[label_indices] = np.random.choice(
            label_names, label_indices.shape[0], p=reassign_probs[:, k]
        )

    adata.obs.loc[adata.obs["is_train"], "labels"] = new_labels_train

    return adata
