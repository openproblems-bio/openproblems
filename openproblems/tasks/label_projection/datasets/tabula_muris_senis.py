from ....data.tabula_muris_senis import load_tabula_muris_senis
from ....tools.decorators import dataset

import numpy as np


@dataset(
    "Tabula Muris Senis Lung (random split)",
    data_url=load_tabula_muris_senis.metadata["data_url"],
    data_reference=load_tabula_muris_senis.metadata["data_reference"],
    dataset_summary="All lung cells from Tabula Muris Senis, a 500k cell-atlas from 18 "
    "organs and tissues across the mouse lifespan. Split into train/test randomly. "
    "Dimensions: 24540 cells, 17985 genes. 39 cell types (avg. 629±999 cells per cell type).",
)
def tabula_muris_senis_lung_random(test=False):
    adata = load_tabula_muris_senis(
        test=test, organ_list=["lung"], method_list=["droplet"]
    )
    adata.obs["labels"] = adata.obs["free_annotation"]
    adata.obs["batch"] = adata.obs["donor_id"]
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
    return adata
