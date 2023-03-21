from os import path
import subprocess
import anndata as ad

name = "pancreas"
output = "dataset.h5ad"
obs_celltype = "celltype"
obs_batch = "tech"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta["executable"],
        "--dataset_id", name,
        "--obs_celltype", obs_celltype,
        "--obs_batch", obs_batch,
        "--layer_counts", "counts",
        "--output", output,
        "--dataset_name", "Pancreas",
        "--data_url", "http://foo.org",
        "--data_reference", "foo2000bar",
        "--dataset_summary", "A short summary.",
        "--dataset_description", "A couple of paragraphs worth of text.",
        "--dataset_organism", "homo_sapiens",
    ],
    check=True
)

print(">> Checking whether file exists", flush=True)
assert path.exists(output)

print(">> Read output anndata", flush=True)
adata = ad.read_h5ad(output)

print(adata)

print(">> Check that output fits expected API", flush=True)
assert adata.X is None
assert "counts" in adata.layers
assert adata.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in adata.obs.columns
if obs_batch:
    assert "batch" in adata.obs.columns

print(">> All tests passed successfully", flush=True)
