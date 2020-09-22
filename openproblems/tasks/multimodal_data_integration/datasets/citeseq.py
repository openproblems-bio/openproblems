from ....data.citeseq import load_citeseq_cmbc
from ....tools.decorators import dataset


@dataset("CITE-seq Cord Blood Mononuclear Cells")
def citeseq_cmbc(test=False):
    return load_citeseq_cmbc(test=test)
