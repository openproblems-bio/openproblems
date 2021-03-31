from ....data.multimodal import scicar
from ....tools.decorators import dataset
from ._helper import _do_dropout


@dataset(
    "sciCAR Mouse Kidney data with evenly distributed dropout in\
    the postive peak counts",
    image="openproblems-python-extras",
)
def scicar_mouse_kidney_dropout(
    test=False, seed=6721, dropout_rate=0.3, cell_fraction=0.8
):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = _do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata
