import subprocess
import scanpy as sc
from os import path

### VIASH START
meta = {
    "resources_dir": "resources_test/label_projection"
}
### VIASH END

input_path = f"{meta['resources_dir']}/pancreas/dataset.h5ad"
output_path = "toy_data.h5ad"

print(">> Running script as test for even")
out = subprocess.check_output([
    meta["executable"],
    "--input", input_path,
    "--output", output_path,
    "--even"
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_path)

print(">> Check that test output fits expected API")
adata = sc.read_h5ad(output_path)
assert (495, 467) == adata.layers["counts"].shape, "processed result data shape {}".format(adata.layers["counts"].shape)

print(">> Runing script as test for specific batch and celltype categories")
out = subprocess.check_output([
    meta["executable"],
    "--input", input_path,
    "--keep_celltype_categories", "acinar:beta",
    "--keep_batch_categories", "celseq:inDrop4:smarter",
    "--output", output_path
]).decode("utf-8")
print(">> Checking whether file exists")
assert path.exists(output_path)

print(">> Check that test output fits expected API")
adata = sc.read_h5ad(output_path)
assert (500, 443) == adata.layers["counts"].shape, "processed result data shape {}".format(adata.layers["counts"].shape)
