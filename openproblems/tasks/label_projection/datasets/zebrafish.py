from ....data.zebrafish import load_zebrafish
from ....tools.decorators import dataset


@dataset("Zebrafish (by labels)")
def zebrafish_labels(test=False):
    adata = load_zebrafish(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["batch"] = adata.obs["lab"]
    adata.obs["is_train"] = adata.obs["lab"] == adata.obs["lab"][0]
    return adata
