import numpy as np
import pyensembl
import os
from ....tools.decorators import method


def _atac_genes_score():
    from pyensembl import EnsemblRelease
    os.system('pyensembl install - -release 86 --species mus_musculus')
    data = EnsemblRelease(86, species='mus_musculus')
    return


@method(
    method_name="BETA",
    paper_name="Target analysis by integration of transcriptome and ChIP-seq data with BETA",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/24263090/",
    paper_year=2013,
    code_version='1.0',
    code_url='',
)
def linear_regression_expotential_decay(adata, decay=1e4):
    adata.obs['atac_rna_cor'] = 1
    return