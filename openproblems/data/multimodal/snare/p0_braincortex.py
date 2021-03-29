from ...utils import loader
from ..utils import subset_joint_data
from .base import load_snare

rna_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5FcDNA%2Ecounts%2Emtx%2Egz"
)
rna_cells_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5FcDNA%2Ebarcodes%2Etsv%2Egz"
)
rna_genes_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5FcDNA%2Egenes%2Etsv%2Egz"
)

atac_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5Fchromatin%2Ecounts%2Emtx%2Egz"
)
atac_cells_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5Fchromatin%2Ebarcodes%2Etsv%2Egz"
)
atac_genes_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074"
    "&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq"
    "%5Fchromatin%2Epeaks%2Etsv%2Egz"
)


@loader
def load_p0_braincortex(test=False):
    """Download sci-CAR mouse kidney data from GEO."""
    if test:
        adata = load_p0_braincortex(test=False)
        adata = subset_joint_data(adata)
        return adata
    return load_snare(
        rna_url, rna_cells_url, rna_genes_url, atac_url, atac_cells_url, atac_genes_url
    )
