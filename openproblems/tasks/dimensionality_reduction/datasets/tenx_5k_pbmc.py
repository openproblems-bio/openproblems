from ....data.tenx import load_tenx_5k_pbmc
from ....tools.decorators import dataset


@dataset(
    "5k Peripheral blood mononuclear cells",
    data_url=load_tenx_5k_pbmc.metadata["data_url"],
    data_reference=load_tenx_5k_pbmc.metadata["data_reference"],
    dataset_summary=(
        "5k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. "
        "Sequenced on 10X v3 chemistry in July 2019 by 10X Genomics."
    ),
)
def tenx_5k_pbmc(test=False):
    return load_tenx_5k_pbmc(test=test)
