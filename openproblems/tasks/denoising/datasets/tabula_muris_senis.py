from ....data.tabula_muris_senis import load_tabula_muris_senis
from ....tools.decorators import dataset
from . import utils


@dataset(
    "Tabula Muris Senis Lung (random split)",
    data_url=load_tabula_muris_senis.metadata["data_url"],
    data_reference=load_tabula_muris_senis.metadata["data_reference"],
    dataset_summary="All lung cells from Tabula Muris Senis, a 500k cell-atlas from 18 "
    "organs and tissues across the mouse lifespan. Split into train/test randomly.",
    image="openproblems-python-extras",
)
def tabula_muris_senis_lung_random(test=False):
    adata = load_tabula_muris_senis(test=test, organ_list=["lung"], method_list=["droplet"])
    adata = utils.split_data(adata)
    return adata
