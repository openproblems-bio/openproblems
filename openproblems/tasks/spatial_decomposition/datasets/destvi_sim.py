from ....tools.decorators import dataset
from ._destvi_utils import generate_synthetic_dataset


@dataset("DestVI simulation (cell)")
def destvi_sim(test=False):

    adata_spatial = generate_synthetic_dataset()
    return adata_spatial
