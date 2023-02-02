from .....data.pancreas import load_pancreas
from .....tools.decorators import dataset
from ..utils import filter_celltypes
from ..utils import precompute_hvg
from typing import Optional


@dataset(
    dataset_name="Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary=(
        "Human pancreatic islet scRNA-seq data from 6 datasets "
        "across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, "
        "and SMARTER-seq)."
    ),
    image="openproblems",
)
def pancreas_batch(
    test: bool = False,
    min_celltype_count: Optional[int] = None,
    n_hvg: Optional[int] = None,
):
    import scanpy as sc

    adata = load_pancreas(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    adata = filter_celltypes(adata, min_celltype_count=min_celltype_count)

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)

    adata.X = adata.layers["log_normalized"]

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni_pca"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni_pca", key_added="uni")

    adata.var_names_make_unique()

    adata.uns["hvg_unint"] = precompute_hvg(adata, n_genes=n_hvg)
    adata.uns["n_genes_pre"] = adata.n_vars
    return adata
