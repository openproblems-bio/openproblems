from .....data.pancreas import load_pancreas
from .....tools.decorators import dataset

import scanpy as sc


@dataset(
    dataset_name="Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreatic islet scRNA-seq data from 6 datasets "
    "across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, "
    "and SMARTER-seq).",
    image="openproblems",
)
def pancreas_batch(test=False):
    adata = load_pancreas(test)
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")

    return adata
