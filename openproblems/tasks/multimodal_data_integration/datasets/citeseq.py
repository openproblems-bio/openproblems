from ....data.multimodal.citeseq import load_citeseq_cbmc
from ....tools.decorators import dataset


@dataset(
    "CITE-seq Cord Blood Mononuclear Cells",
    data_url=load_citeseq_cbmc.metadata["data_url"],
    data_reference=load_citeseq_cbmc.metadata["data_reference"],
    dataset_summary="8k cord blood mononuclear cells sequenced by CITEseq, a multimodal"
    " addition to the 10x scRNA-seq platform that allows simultaneous measurement of "
    "RNA and protein.",
)
def citeseq_cbmc(test=False):
    return load_citeseq_cbmc(test=test)
