from ....data.multimodal import scicar
from ....tools.decorators import dataset


@dataset(
    "sciCAR Mouse Kidney with cell clusters",
    data_url=scicar.load_scicar_mouse_kidney.metadata["data_url"],
    data_reference=scicar.load_scicar_mouse_kidney.metadata["data_reference"],
    dataset_summary="11k cells from adult mouse kidney sequenced "
    "by sci-CAR, a combinatorial indexing-based co-assay that jointly profiles "
    "chromatin accessibility and mRNA.",
)
def scicar_mouse_kidney(test=False):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata
