from ....data.tnbc_wu2021 import load_tnbc_data
from ....tools.decorators import dataset


@dataset(
    "~43k Cells from Triple Negative Breast Cancer"
    "Wu, S.Z., Al-Eryani, G., Roden, D.L., Junankar, S., Harvey, K., Andersson, A., Thennavan, A., Wang, C., Torpy, J.R., Bartonicek, N. and Wang, T., 2021. A single-cell and spatially resolved atlas of human breast cancers. Nature genetics, 53(9), pp.1334-1347.",
    data_url=load_tnbc_data.metadata["data_url"],
    dataset_summary="TODO",
)
def tnbc_data(test=False):
    return load_tnbc_data(test=test)
