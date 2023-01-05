from typing import Optional


def filter_celltypes(
    adata, test: bool = False, min_celltype_count: Optional[int] = None
):

    if test:
        min_celltype_count = min_celltype_count or 2
    else:
        min_celltype_count = min_celltype_count or 50

    celltype_counts = adata.obs["labels"].value_counts()
    keep_celltypes = celltype_counts[celltype_counts >= min_celltype_count].index
    keep_cells = adata.obs["labels"].isin(keep_celltypes)
    return adata[keep_cells].copy()
