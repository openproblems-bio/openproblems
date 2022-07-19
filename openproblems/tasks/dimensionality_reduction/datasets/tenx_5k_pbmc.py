from ....data.tenx_5k_pbmc import load_tenx_5k_pbmc
from ....tools.decorators import dataset


@dataset(
    "5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor. "
    "10x Genomics; July 24, 2019.",
    data_url=load_tenx_5k_pbmc.metadata["data_url"],
    data_reference=load_tenx_5k_pbmc.metadata["data_reference"],
    dataset_summary="TODO",
)
def tenx_5k_pbmc(test=False):
    return load_tenx_5k_pbmc(test=test)
