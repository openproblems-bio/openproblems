from ....data.multimodal import snare
from ....tools.decorators import dataset
from ._helper import _do_dropout


@dataset(
    "SNARE-seq P0 brain cortex data with evenly distributed dropout in\
    the postive peak counts",
    image="openproblems-python-extras",
)
def snare_p0_braincortex_dropout(
    test=False, seed=243287, dropout_rate=0.3, cell_fraction=0.8
):
    adata = snare.load_p0_braincortex(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = _do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata
