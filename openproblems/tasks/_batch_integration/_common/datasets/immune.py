from .....data.immune_cells import load_immune
from .....tools.decorators import dataset
from ..utils import filter_celltypes
from typing import Optional


@dataset(
    dataset_name="Immune (by batch)",
    data_url=load_immune.metadata["data_url"],
    data_reference=load_immune.metadata["data_reference"],
    dataset_summary="Human immune cells from peripheral blood and bone marrow "
    "taken from 5 datasets comprising 10 batches across technologies (10X, "
    "Smart-seq2).",
    image="openproblems",
)
def immune_batch(test: bool = False, min_celltype_count: Optional[int] = None):
    import scanpy as sc

    adata = load_immune(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["final_annotation"]

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
