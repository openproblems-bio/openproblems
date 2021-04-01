from ....data.multimodal import share
from ....tools.decorators import dataset
from .utils import do_dropout


@dataset(
    "SNARE-seq mouse skin data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_mouse_skin_dropout(test=False, seed=683, dropout_rate=0.3, cell_fraction=0.8):
    adata = share.load_share_mouse_skin(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "SNARE-seq mouse brain data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_mouse_brain_dropout(
    test=False, seed=28479, dropout_rate=0.3, cell_fraction=0.8
):
    adata = share.load_share_mouse_brain(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "SNARE-seq mouse lung data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_mouse_lung_dropout(
    test=False, seed=17973, dropout_rate=0.3, cell_fraction=0.8
):
    adata = share.load_share_mouse_lung(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "SNARE-seq GM12878 rep1 data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_gm12878_rep1_dropout(
    test=False, seed=346332, dropout_rate=0.3, cell_fraction=0.8
):
    adata = share.load_share_gm12878_rep1(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "SNARE-seq GM12878 rep2 data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_gm12878_rep2_dropout(
    test=False, seed=34637, dropout_rate=0.3, cell_fraction=0.8
):
    adata = share.load_share_gm12878_rep2(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata


@dataset(
    "SNARE-seq GM12878 rep3 data with more peak count dropouts",
    image="openproblems-python-extras",
)
def share_gm12878_rep3_dropout(
    test=False, seed=368342, dropout_rate=0.3, cell_fraction=0.8
):
    adata = share.load_share_gm12878_rep3(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    adata = do_dropout(
        adata, seed, dropout_rate=dropout_rate, cell_fraction=cell_fraction
    )
    return adata
