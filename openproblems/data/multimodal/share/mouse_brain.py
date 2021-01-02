from ...utils import loader
from ..utils import subset_joint_data
from .base import load_share

rna_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156610"
    "&format=file&file=GSM4156610%5Fbrain%2Erna%2Ecounts%2Etxt%2Egz"
)

atac_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156599"
    "&format=file&file=GSM4156599%5Fbrain%2Ecounts%2Etxt%2Egz"
)
atac_cells_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156599"
    "&format=file&file=GSM4156599%5Fbrain%2Ebarcodes%2Etxt%2Egz"
)
atac_genes_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156599"
    "&format=file&file=GSM4156599%5Fbrain%2Epeaks%2Ebed%2Egz"
)


@loader
def load_share_mouse_brain(test=False):
    """Download sci-CAR mouse kidney data from GEO."""
    if test:
        adata = load_share_mouse_brain(test=False)
        adata = subset_joint_data(adata)
        return adata
    return load_share(rna_url, atac_url, atac_cells_url, atac_genes_url)
