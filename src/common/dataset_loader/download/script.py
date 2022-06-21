print("Importing libraries")
import scanpy as sc
import tempfile
import os
import urllib

_FAKE_HEADERS = [("User-Agent", "Mozilla/5.0")]

## VIASH START
par = {
    "url": "https://ndownloader.figshare.com/files/24539828",
    "name": "pancreas",
    "obs_celltype": "celltype",
    "obs_batch": "tech",
    "output": "test_data.h5ad"
}
## VIASH END

with tempfile.TemporaryDirectory() as tempdir:
    print("Downloading", par['url'], flush=True)
    filepath = os.path.join(tempdir, "dataset.h5ad")

    with open(filepath, "wb") as filehandle:
        opener = urllib.request.build_opener()
        opener.addheaders = _FAKE_HEADERS
        urllib.request.install_opener(opener)
        with urllib.request.urlopen(par["url"]) as urlhandle:
            filehandle.write(urlhandle.read())
    # scprep.io.download.download_url(par['url'], filepath)

    print("Reading file")
    adata = sc.read_h5ad(filepath)

print("Copying .layers['counts'] to .X")
if "counts" in adata.layers:
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]

print("Setting .uns['dataset_id']")
adata.uns["dataset_id"] = par["name"]

print("Setting .obs['celltype']")
if par["obs_celltype"]:
    if par["obs_celltype"] in adata.obs:
        adata.obs["celltype"] = adata.obs[par["obs_celltype"]]
    else:
        print(f"Warning: key '{par['obs_celltype']}' could not be found in adata.obs.")

print("Setting .obs['batch']")
if par["obs_batch"]:
    if par["obs_batch"] in adata.obs:
        adata.obs["batch"] = adata.obs[par["obs_batch"]]
    else:
        print(f"Warning: key '{par['obs_batch']}' could not be found in adata.obs.")

print("Remove cells or genes with 0 counts")
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=2)

print("Writing adata to file")
adata.write(par["output"], compression="gzip")
