## VIASH START
par = {
    "url": "https://ndownloader.figshare.com/files/24539828",
    "output": "/tmp/output.h5ad"
}
## VIASH END
print("Import libraries")
import os
import scanpy as sc
import scprep
import tempfile


with tempfile.TemporaryDirectory() as tempdir:
    filepath = os.path.join(tempdir, "pancreas.h5ad")
    print("Download data from '--url'")
    scprep.io.download.download_url(par['url'], filepath)
    adata = sc.read(filepath)
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]

print("Writing data in {output}".format(output=par['output']))
adata.write(par['output'])
