from ....data.zebrafish import load_zebrafish
from ....tools.decorators import dataset


@dataset(
    "Zebrafish",
    data_url=load_zebrafish.metadata["data_url"],
    data_reference=load_zebrafish.metadata["data_reference"],
    dataset_summary="90k cells from zebrafish embryos throughout the first day of "
    "development, with and without a knockout of chordin, an important developmental "
    "gene. Dimensions: 26022 cells, 25258 genes. 24 cell types "
    "(avg. 1084±1156 cells per cell type).",
)
def zebrafish_labs(test=False):
    adata = load_zebrafish(test=test)
    return adata
