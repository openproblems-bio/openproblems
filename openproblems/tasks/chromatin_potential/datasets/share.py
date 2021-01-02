from ....data.multimodal import share
from ....tools.decorators import dataset


@dataset("SHARE-seq mouse skin data with cell clusters")
def share_mouse_skin(test=False):
    adata = share.load_share_mouse_skin(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata


@dataset("SHARE-seq mouse brain data with cell clusters")
def share_mouse_brain(test=False):
    adata = share.load_share_mouse_brain(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata


@dataset("SHARE-seq mouse lung data with cell clusters")
def share_mouse_lung(test=False):
    adata = share.load_share_mouse_lung(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata


@dataset("SHARE-seq GM12878 rep1 data with cell clusters")
def share_gm12878_rep1(test=False):
    adata = share.load_share_gm12878_rep1(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    return adata


@dataset("SHARE-seq GM12878 rep2 data with cell clusters")
def share_gm12878_rep2(test=False):
    adata = share.load_share_gm12878_rep2(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    return adata


@dataset("SHARE-seq GM12878 rep3 data with cell clusters")
def share_gm12878_rep3(test=False):
    adata = share.load_share_gm12878_rep3(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "97"
    return adata
