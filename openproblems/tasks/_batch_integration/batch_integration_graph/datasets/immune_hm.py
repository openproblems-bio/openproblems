from .....data.immune_cells_human_mouse import load_immune_hm
from .....tools.decorators import dataset


@dataset(
    dataset_name="Immune (human/mouse)(by batch)",
    data_url=load_immune_hm.metadata["data_url"],
    data_reference=load_immune_hm.metadata["data_reference"],
    dataset_summary="Human and Mouse immune cells from peripheral blood and bone marrow "
    "taken from multiple datasets across technologies (10X, "
    "Smart-seq2).",
    image="openproblems",
)
def immune_hm_batch(test=False):
    import scanpy as sc

    adata = load_immune_hm(test)
    adata.uns["organism"] = "human"
    adata.obs["labels"] = adata.obs["final_annotation"]

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
