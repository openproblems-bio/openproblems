from ....data.multimodal import scicar
from ....tools.decorators import dataset
from .utils import do_dropout


@dataset(
    "sciCAR mouse kidney data with more peak count dropouts",
    image="openproblems-python-extras",
)
def scicar_mouse_kidney_dropout(
    test=False, seed=6721, dropout_rate=0.3, cell_fraction=0.8
):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "sciCAR mouse kidney data with only 90% reads",
    image="openproblems-python-extras",
)
def scicar_mouse_kidney_split(
    test=False, seed=6721, dropout_rate=0.3, cell_fraction=0.8
):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata
