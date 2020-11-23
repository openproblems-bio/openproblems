  
from ....data._10x_5k_pbmc import load_10x_5k_pbmc
from ....tools.decorators import dataset


@dataset("10X PBMCs v3 Chemistry")
def _10x_5k_pbmc(test=False):
    return load_10x_5k_pbmc(test=test)
