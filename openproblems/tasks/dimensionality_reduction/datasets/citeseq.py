from ....data.multimodal.citeseq import load_citeseq_cbmc
from ....tools.decorators import dataset
from .preprocessing import preprocess_scanpy


@dataset("CITE-seq Cord Blood Mononuclear Cells")
def citeseq_cbmc(test=False):

    adata = load_citeseq_cbmc(test=test)

    return preprocess_scanpy(adata)
