from . import utils
from ....data.pbmc import load_pbmc
from ....tools.decorators import dataset


@dataset("1k PBMCs from a Healthy Donor (10x/v3)")
def pbmc(test=False):
    adata = load_pbmc(test=test)
    adata = utils.split_data(adata)
    return adata
