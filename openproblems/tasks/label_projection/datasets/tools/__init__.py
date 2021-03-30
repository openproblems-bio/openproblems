import numpy as np
import pandas as pd


def add_label_noise(adata, noise_prob):
    """Inject random label noise in the dataset by permuting a fraction
    of the labels in the training set.

    By adding different levels of label noise metrics can be evaluated to show
    generalization trends from training data even if ground truth is uncertain.

    Params
    -------
    adata: AnnData, a dataset with the required fields for the label_projection task.

    noise_prob: Float, the probability of label noise in the training data.

    Returns
    -------
    AnnData: dataset where training labels have been permuted by specified probability.
    """

    old_labels = adata.obs["labels"].pipe(np.array)
    new_labels = adata.obs["labels"].pipe(np.array)

    label_names = np.unique(new_labels)

    n_labels = label_names.shape[0]

    reassign_probs = (noise_prob / (n_labels - 1)) * np.ones((n_labels, n_labels))

    np.fill_diagonal(reassign_probs, 1 - noise_prob)

    for k, label in enumerate(label_names):
        label_indices = np.where(old_labels == label)[0]
        new_labels[label_indices] = np.random.choice(
            label_names, label_indices.shape[0], p=reassign_probs[:, k]
        )
