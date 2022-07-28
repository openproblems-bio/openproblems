from ....data.tabula_muris_senis import load_tabula_muris_senis
from ....tools.decorators import dataset

import numpy as np


@dataset(
    "Tabula Muris Senis Lung (random split)",
    data_url=load_tabula_muris_senis.metadata["data_url"],
    data_reference=load_tabula_muris_senis.metadata["data_reference"],
    dataset_summary="TODO",
)
def tabula_muris_senis_lung_random(test=False):
    adata = load_tabula_muris_senis(test=test, organ_list=["lung"])
    adata.obs["labels"] = adata.obs["free_annotation"]
    adata.obs["batch"] = adata.obs["mouse.id"]
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
    return adata
