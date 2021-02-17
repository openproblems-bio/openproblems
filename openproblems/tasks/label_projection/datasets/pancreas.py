import numpy as np

from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset


@dataset("Pancreas (by batch)")
def pancreas_batch(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]

    # Assign training/test
    test_batches = adata.obs["tech"].dtype.categories[[-3, -1]]
    adata.obs["is_train"] = [
        False if adata.obs["tech"][idx] in test_batches else True
        for idx in adata.obs_names
    ]

    return adata


@dataset("Pancreas (random split)")
def pancreas_random(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]

    # Assign training/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    return adata
