from ....data.utils import filter_genes_cells
from ....tools.decorators import dataset
from ._destvi_utils import generate_synthetic_dataset_destvi


@dataset(
    "Destvi",
    data_url="https://www.nature.com/articles/s41587-022-01272-8",
    dataset_summary="Synthetic dataset from destvi manuscript.",
)
def destvi_dataset():
    merged_adata = generate_synthetic_dataset_destvi()
    filter_genes_cells(merged_adata)
    return merged_adata
