from os import path
import subprocess
import anndata as ad

input_id = "pancreas"
dataset_id = "openproblems_v1/" + input_id
output = "dataset.h5ad"
obs_cell_type = "celltype"
obs_batch = "tech"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta["executable"],
        "--input_id", input_id,
        "--dataset_id", dataset_id,
        "--obs_cell_type", obs_cell_type,
        "--obs_batch", obs_batch,
        "--layer_counts", "counts",
        "--output", output,
        "--dataset_name", "Pancreas",
        "--dataset_url", "http://foo.org",
        "--dataset_reference", "foo2000bar",
        "--dataset_summary", "A short summary.",
        "--dataset_description", "A couple of paragraphs worth of text.",
        "--dataset_organism", "homo_sapiens",
    ],
    stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.")
    exit(out.returncode)

print(">> Checking whether file exists", flush=True)
assert path.exists(output), "Output does not exist"

print(">> Read output anndata", flush=True)
adata = ad.read_h5ad(output)

print(adata)

print(">> Check that output fits expected API", flush=True)
assert adata.X is None, "adata.X should be None/empty"
assert "counts" in adata.layers, "Counts layer not found in output layers"
assert adata.uns["dataset_id"] == dataset_id, f"Expected {dataset_id} as value"
if obs_cell_type:
    assert "cell_type" in adata.obs.columns, "'cell_type' column not found in obs of anndata output"
if obs_batch:
    assert "batch" in adata.obs.columns, "'batch' column not found in obs of anndata output"

print(">> All tests passed successfully", flush=True)
