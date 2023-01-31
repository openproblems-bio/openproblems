from .....data.lung import load_lung
from .....tools.decorators import dataset


@dataset(
    dataset_name="Lung (Viera Braga et al.)",
    data_url=load_lung.metadata["data_url"],
    data_reference=load_lung.metadata["data_reference"],
    dataset_summary="Human lung scRNA-seq data from 3 datasets with 32,472 cells."
    "From Vieira Braga et al. Technologies: 10X and Drop-seq.",
    image="openproblems",
)
def lung_batch(test=False):
    import scanpy as sc

    adata = load_lung(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["cell_type"]
    # No need to rename batch column as it already exists

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
