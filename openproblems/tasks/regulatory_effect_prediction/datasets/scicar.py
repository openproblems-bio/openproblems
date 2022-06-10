from ....data.multimodal import scicar
from ....tools.decorators import dataset


@dataset(
    "sciCAR Mouse Kidney with cell clusters",
    data_url=scicar.load_scicar_mouse_kidney.metadata["data_url"],
)
def scicar_mouse_kidney(test=False):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata
