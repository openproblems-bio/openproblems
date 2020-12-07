import scanpy as sc
from .utils import loader

url = "https://ndownloader.figshare.com/files/24130931"


@loader
def load_senis_muris(test=False):
    """Download Tabula Senis."""
    adata = sc.read(url=url)
    if test:
        sc.pp.subsample(adata, fraction=0.1)
    return adata
