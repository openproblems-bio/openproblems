from ....data.pbmc import load_pbmc
from ....tools.decorators import dataset
from . import utils


@dataset(
    "1k PBMCs from a Healthy Donor (10x/v3)",
    data_url=load_pbmc.metadata["data_url"],
    data_reference=load_pbmc.metadata["data_reference"],
    dataset_summary="TODO",
    image="openproblems-python-extras",
)
def pbmc(test=False):
    adata = load_pbmc(test=test)
    adata = utils.split_data(adata)
    return adata
