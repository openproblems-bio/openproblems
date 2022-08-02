from ....data.tenx import load_tenx_1k_pbmc
from ....tools.decorators import dataset
from . import utils


@dataset(
    "1k Peripheral blood mononuclear cells",
    data_url=load_tenx_1k_pbmc.metadata["data_url"],
    data_reference=load_tenx_1k_pbmc.metadata["data_reference"],
    dataset_summary=(
        "1k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. "
        "Sequenced on 10X v3 chemistry in November 2018 by 10X Genomics."
    ),
    image="openproblems-python-extras",
)
def pbmc(test=False):
    adata = load_tenx_1k_pbmc(test=test)
    adata = utils.split_data(adata)
    return adata
