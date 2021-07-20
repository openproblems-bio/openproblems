## VIASH START
par = {
    "url": "https://ndownloader.figshare.com/files/24539828",  # pancreas data
    "name": "pancreas",
    "output": "./src/common/data_loader/resources/pancreas.h5ad",
    "lognorm_available": "adata.X"
}
## VIASH END

print("Importing libraries")
import sys
import scanpy as sc
import tempfile
import os
import scprep


with tempfile.TemporaryDirectory() as tempdir:
    URL = par['url']
    print("Downloading", URL)
    sys.stdout.flush()
    filepath = os.path.join(tempdir, "pancreas.h5ad")
    scprep.io.download.download_url(URL, filepath)

    print("Read file")
    adata = sc.read(filepath)
    adata.uns["name"] = par["name"]

    # Rearrange preprocessing
    if par['lognorm_available'] == 'adata.X':
        # save lognorm counts
        adata.layers['logcounts'] = adata.X

    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        del adata.layers["counts"]

    # Ensure there are no cells or genes with 0 counts
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)


print("Writing adata to file")
adata.write(par["output"], compression="gzip")
