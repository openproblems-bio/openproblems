from typing import Optional
from scanpy.pp import highly_variable_genes


def filter_celltypes(adata, min_celltype_count: Optional[int] = None):

    min_celltype_count = min_celltype_count or 50

    celltype_counts = adata.obs["labels"].value_counts()
    keep_celltypes = celltype_counts[celltype_counts >= min_celltype_count].index
    keep_cells = adata.obs["labels"].isin(keep_celltypes)
    return adata[keep_cells].copy()


def precomp_hvg(adata, n_genes):
    hvg_unint = highly_variable_genes(
        adata,
        n_top_genes=n_genes,
        layer="log_normalized",
        flavor="cell_ranger",
        batch_key="batch",
        inplace=False,
    )
    return list(hvg_unint[hvg_unint.highly_variable].index)
