import numpy as np
from ....data.pancreas import load_pancreas


def pancreas_batch(test=False):
    adata = load_pancreas(test=test)

    # Assign training/test
    test_batches = ["inDrop4", "smartseq2"]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]

    return adata


def pancreas_random(test=False):
    adata = load_pancreas(test=test)

    # Assign training/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    return adata
