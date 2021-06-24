## VIASH START
par = {
    "url": "https://ndownloader.figshare.com/files/24974582",  # PBMC data
    "output": "test_data.h5ad"
}
resources_dir = "./utils/"
## VIASH END

print("Importing libraries")
import scanpy as sc
import tempfile
import os
import scprep

# # adding resources dir to system path
# # the resources dir contains all files listed in the '.functionality.resources' part of the
# # viash config, amongst which is the 'utils.py' file we need.
import sys

sys.path.append(resources_dir)
# # importing helper functions from common utils.py file in resources dir
from utils import filter_genes_cells

with tempfile.TemporaryDirectory() as tempdir:
    URL = par['url']
    print("Downloading", URL)
    sys.stdout.flush()
    filepath = os.path.join(tempdir, "pancreas.h5ad")
    scprep.io.download.download_url(URL, filepath)

    print("Read file")
    adata = sc.read(filepath)

    # Remove preprocessing
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        del adata.layers["counts"]

    # Ensure there are no cells or genes with 0 counts
    filter_genes_cells(adata)


print("Writing adata to file")
adata.write(par["output"], compression="gzip")
