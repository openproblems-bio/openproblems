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
    "layer_counts": "counts",
    "output": "test_data.h5ad",
    "layer_counts_output": "counts"
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

    print("Reading file")
    adata = sc.read_h5ad(filepath)

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

print("Setting .obs['tissue']")
if par["obs_tissue"]:
    if par["obs_tissue"] in adata.obs:
        adata.obs["tissue"] = adata.obs[par["obs_tissue"]]
    else:
        print(f"Warning: key '{par['obs_tissue']}' could not be found in adata.obs.")

print("Remove cells or genes with 0 counts")
if par["layer_counts"] and par["layer_counts"] in adata.layers:
    print(f"  Temporarily copying .layers['{par['layer_counts']}'] to .X")
    adata.X = adata.layers[par["layer_counts"]]
    del adata.layers[par["layer_counts"]]

print("  Removing empty genes")
sc.pp.filter_genes(adata, min_cells=1)
print("  Removing empty cells")
sc.pp.filter_cells(adata, min_counts=2)

if par["layer_counts_output"]:
    print(f"  Copying .X back to .layers['{par['layer_counts_output']}']")
    adata.layers[par["layer_counts_output"]] = adata.X
    del adata.X

print("Writing adata to file")
adata.write_h5ad(par["output"], compression="gzip")
