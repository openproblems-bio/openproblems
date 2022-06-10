from ....data.multimodal.citeseq import load_citeseq_cbmc
from ....tools.decorators import dataset


@dataset(
    "CITE-seq Cord Blood Mononuclear Cells",
    data_url=load_citeseq_cbmc.metadata["data_url"],
)
def citeseq_cbmc(test=False):
    return load_citeseq_cbmc(test=test)
