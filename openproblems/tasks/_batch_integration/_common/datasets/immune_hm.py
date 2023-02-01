from .....data.immune_cells_human_mouse import load_immune_human_mouse
from .....tools.decorators import dataset
from ..utils import filter_celltypes
from typing import Optional


@dataset(
    dataset_name="Immune (human/mouse)",
    data_url=load_immune_human_mouse.metadata["data_url"],
    data_reference=load_immune_human_mouse.metadata["data_reference"],
    dataset_summary="Human and Mouse immune cells from peripheral blood and bone marrow"
    "taken from 8 datasets with 23 batches across 5 technologies (10X 3 v2,"
    "10x 3 v3, microwell seq and Smart-seq2). The data contains 97861 cells in total."
    "Mouse gene names are mapped to human gene names by capitalisation.",
    image="openproblems",
)
def immune_human_mouse_batch(
    test: bool = False, min_celltype_count: Optional[int] = None
):
    import scanpy as sc

    adata = load_immune_human_mouse(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["final_annotation"]
    # No need to rename batch column as it already exists

    adata = filter_celltypes(adata, min_celltype_count=min_celltype_count)

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    adata.X = adata.layers["log_normalized"]

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni_pca"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni_pca", key_added="uni")
    adata.var_names_make_unique()
    return adata
