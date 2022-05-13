## VIASH START
par = {
    "url": "https://ndownloader.figshare.com/files/24539828",
    "output": "/tmp/output.h5ad"
}
## VIASH END
import os
from os.path import exists
import scanpy as sc
import scprep
import tempfile


def load_data(url):
    """Download pancreas data from url."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "pancreas.h5ad")
        scprep.io.download.download_url(url, filepath)
        adata = sc.read(filepath)
        adata.X = adata.layers["counts"]
        del adata.layers["counts"]
    return adata


def save_data(adata, output):
    adata.write(output)
    return None


def run_script(par):
    adata = load_data(par['url'])
    save_data(adata, par['output'])
    return None


if ("par" in locals()):
    run_script(par)
