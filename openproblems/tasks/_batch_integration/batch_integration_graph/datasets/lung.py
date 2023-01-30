from .....data.lung import load_lung
from .....tools.decorators import dataset


@dataset(
    dataset_name="Lung (by batch)",
    data_url=load_lung.metadata["data_url"],
    data_reference=load_lung.metadata["data_reference"],
    dataset_summary="Human pancreatic islet scRNA-seq data from 6 datasets "
    "across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, "
    "and SMARTER-seq).",
    image="openproblems",
)
def lung_batch(test=False):
    import scanpy as sc

    adata = load_lung(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["cell_type"]

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
    return adata
