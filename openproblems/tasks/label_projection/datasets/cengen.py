from ....data.cengen import load_cengen
from ....tools.decorators import dataset

import numpy as np


@dataset(
    "CeNGEN (by batch)",
    data_url=load_cengen.metadata["data_url"],
    dataset_summary="TODO",
)
def cengen_batch(test=False):
    adata = load_cengen(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["batch"] = adata.obs["experiment_code"]

    # Assign training/test
    test_batches = [adata.obs.batch.value_counts().index[-1]]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]

    return adata


@dataset(
    "CeNGEN (random split)",
    data_url=load_cengen.metadata["data_url"],
    dataset_summary="TODO",
)
def cengen_random(test=False):
    adata = load_cengen(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["batch"] = adata.obs["experiment_code"]

    # Assign training/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    return adata
