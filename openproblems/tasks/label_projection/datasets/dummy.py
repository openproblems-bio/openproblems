import numpy as np
from ....data.dummy import load_dummy


def dummy(test=False):
    adata = load_dummy()
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True
    )
    return adata
