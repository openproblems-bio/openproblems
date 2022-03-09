from ...utils import loader
from ..utils import subset_joint_data
from .base import load_share

rna_url1 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156601"
    "&format=file&file=GSM4156601%5FGM12878%2Erep1%2Erna%2Ecounts%2Etxt%2Egz"
)

atac_url1 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156590"
    "&format=file&file=GSM4156590%5FGM12878%2Erep1%2Ecounts%2Etxt%2Egz"
)
atac_cells_url1 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156590"
    "&format=file&file=GSM4156590%5FGM12878%2Erep1%2Ebarcodes%2Etxt%2Egz"
)
atac_genes_url1 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156590"
    "&format=file&file=GSM4156590%5FGM12878%2Erep1%2Epeaks%2Ebed%2Egz"
)

rna_url2 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156602"
    "&format=file&file=GSM4156602%5FGM12878%2Erep2%2Erna%2Ecounts%2Etxt%2Egz"
)

atac_url2 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156591"
    "&format=file&file=GSM4156591%5FGM12878%2Erep2%2Ecounts%2Etxt%2Egz"
)
atac_cells_url2 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156591"
    "&format=file&file=GSM4156591%5FGM12878%2Erep2%2Ebarcodes%2Etxt%2Egz"
)
atac_genes_url2 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156591"
    "&format=file&file=GSM4156591%5FGM12878%2Erep2%2Epeaks%2Ebed%2Egz"
)

rna_url3 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156603"
    "&format=file&file=GSM4156603%5FGM12878%2Erep3%2Erna%2Ecounts%2Etxt%2Egz"
)

atac_url3 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156592"
    "&format=file&file=GSM4156592%5FGM12878%2Erep3%2Ecounts%2Etxt%2Egz"
)
atac_cells_url3 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156592"
    "&format=file&file=GSM4156592%5FGM12878%2Erep3%2Ebarcodes%2Etxt%2Egz"
)
atac_genes_url3 = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4156592"
    "&format=file&file=GSM4156592%5FGM12878%2Erep3%2Epeaks%2Ebed%2Egz"
)


@loader
def load_share_gm12878_rep1(test=False):
    """Download sci-CAR mouse kidney data from GEO."""
    if test:
        adata = load_share_gm12878_rep1(test=False)
        adata = subset_joint_data(adata)
        return adata
    return load_share(rna_url1, atac_url1, atac_cells_url1, atac_genes_url1)


@loader
def load_share_gm12878_rep2(test=False):
    """Download sci-CAR mouse kidney data from GEO."""
    if test:
        adata = load_share_gm12878_rep2(test=False)
        adata = subset_joint_data(adata)
        return adata
    return load_share(rna_url2, atac_url2, atac_cells_url2, atac_genes_url2)


@loader
def load_share_gm12878_rep3(test=False):
    """Download sci-CAR mouse kidney data from GEO."""
    if test:
        adata = load_share_gm12878_rep3(test=False)
        adata = subset_joint_data(adata)
        return adata
    return load_share(rna_url3, atac_url3, atac_cells_url3, atac_genes_url3)
