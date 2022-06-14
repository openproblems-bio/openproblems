from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset
from .tools import add_label_noise

import numpy as np


@dataset(
    "Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_batch(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign training/test
    test_batches = adata.obs["batch"].dtype.categories[[-3, -1]]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]

    return adata


@dataset(
    "Pancreas (random split)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_random(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign training/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    return adata


@dataset(
    r"Pancreas (random split with 20% label noise",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
)
def pancreas_random_label_noise(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign trainin/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    # Inject label noise
    adata = add_label_noise(adata, noise_prob=0.2)

    return adata
