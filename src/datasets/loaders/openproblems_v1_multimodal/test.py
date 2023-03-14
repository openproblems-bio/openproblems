from os import path
import subprocess
import anndata as ad

name = "scicar_mouse_kidney"
obs_celltype = "cell_name"
obs_batch = "replicate"
obs_tissue = None

output_mod1_file = "output_mod1.h5ad"
output_mod2_file = "output_mod2.h5ad"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta["executable"],
        "--dataset_id", name,
        "--obs_celltype", obs_celltype,
        "--obs_batch", obs_batch,
        "--layer_counts", "counts",
        "--output_mod1", output_mod1_file,
        "--output_mod2", output_mod2_file,
        "--dataset_name", "Pancreas",
        "--data_url", "http://foo.org",
        "--data_reference", "foo2000bar",
        "--dataset_summary", "A short summary.",
        "--dataset_description", "A couple of paragraphs worth of text.",
        "--dataset_organism", "homo_sapiens",
    ],
    check=True
)

print(">> Checking whether files exist", flush=True)
assert path.exists(output_mod1_file)
assert path.exists(output_mod2_file)

print(">> Read output anndata", flush=True)
output_mod1 = ad.read_h5ad(output_mod1_file)
output_mod2 = ad.read_h5ad(output_mod2_file)

print(f"output_mod1: {output_mod1}", flush=True)
print(f"output_mod2: {output_mod2}", flush=True)

print(">> Check that output mod1 fits expected API", flush=True)
assert output_mod1.X is None
assert "counts" in output_mod1.layers
assert output_mod1.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in output_mod1.obs.columns
if obs_batch:
    assert "batch" in output_mod1.obs.columns
if obs_tissue:
    assert "tissue" in output_mod1.obs.columns

print(">> Check that output mod2 fits expected API", flush=True)
assert output_mod2.X is None
assert "counts" in output_mod2.layers
assert output_mod2.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in output_mod2.obs.columns
if obs_batch:
    assert "batch" in output_mod2.obs.columns
if obs_tissue:
    assert "tissue" in output_mod2.obs.columns

print(">> All tests passed successfully", flush=True)
