from ....data.multimodal import snare
from ....tools.decorators import dataset


@dataset("SNARE-seq P0 brain cortex")
def snare_p0_braincortex(test=False):
    adata = snare.load_p0_braincortex(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"
    return adata
