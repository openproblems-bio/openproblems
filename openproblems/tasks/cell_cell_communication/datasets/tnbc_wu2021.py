from ....data.tnbc_wu2021 import load_tnbc_data
from ....tools.decorators import dataset


@dataset(
    "~43k Cells from Triple Negative Breast Cancer"
    "Wu et al., 2021. Nature genetics, 53(9), pp.1334-1347.",
    data_url=load_tnbc_data.metadata["data_url"],
    dataset_summary="TODO",
)
def tnbc_data(test=False):
    return load_tnbc_data(test=test)
