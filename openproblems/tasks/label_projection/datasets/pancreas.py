import numpy as np
from .....data.pancreas import load_pancreas


def pancreas_batch(test=False):
    adata = load_pancreas(test=test)

    # Assign training/test
    test_batches = ["inDrop4", "smartseq2"]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]

    return adata
