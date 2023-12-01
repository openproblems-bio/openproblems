from os import path
import subprocess
import anndata as ad

name = "scicar_mouse_kidney"
obs_cell_type = "cell_name"
obs_batch = "replicate"
obs_tissue = None

output_mod1_file = "output_mod1.h5ad"
output_mod2_file = "output_mod2.h5ad"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta["executable"],
        "--dataset_id", name,
        "--obs_cell_type", obs_cell_type,
        "--obs_batch", obs_batch,
        "--layer_counts", "counts",
        "--output_mod1", output_mod1_file,
        "--output_mod2", output_mod2_file,
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

print(">> Checking whether files exist", flush=True)
assert path.exists(output_mod1_file), "Output mod1 file does not exist"
assert path.exists(output_mod2_file), "Output mod2 file does not exist"

print(">> Read output anndata", flush=True)
output_mod1 = ad.read_h5ad(output_mod1_file)
output_mod2 = ad.read_h5ad(output_mod2_file)

print(f"output_mod1: {output_mod1}", flush=True)
print(f"output_mod2: {output_mod2}", flush=True)

print(">> Check that output mod1 fits expected API", flush=True)
assert output_mod1.X is None, ".X is not None/empty in mod 1 output"
assert "counts" in output_mod1.layers, "'counts' not found in mod 1 output layers" 
assert output_mod1.uns["dataset_id"] == name, f"Expected: {name} as value for dataset_id in mod 1 output uns"
if obs_cell_type:
    assert "cell_type" in output_mod1.obs.columns, "cell_type column not found in mod 1 output obs"
if obs_batch:
    assert "batch" in output_mod1.obs.columns, "batch column not found in mod 1 output obs"
if obs_tissue:
    assert "tissue" in output_mod1.obs.columns, "tissue column not found in mod 1 output obs"

print(">> Check that output mod2 fits expected API", flush=True)
assert output_mod2.X is None, ".X is not None/empty in mod 2 output"
assert "counts" in output_mod2.layers, "'counts' not found in mod 2 output layers"
assert output_mod2.uns["dataset_id"] == name, f"Expected: {name} as value for dataset_id in mod 2 output uns"
if obs_cell_type:
    assert "cell_type" in output_mod2.obs.columns, "cell_type column not found in mod 2 output obs"
if obs_batch:
    assert "batch" in output_mod2.obs.columns, "batch column not found in mod 2 output obs"
if obs_tissue:
    assert "tissue" in output_mod2.obs.columns, "tissue column not found in mod 2 output obs"

print(">> All tests passed successfully", flush=True)
