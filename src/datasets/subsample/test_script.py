import subprocess
import scanpy as sc
from os import path

### VIASH START
meta = {
    "resources_dir": "resources_test/common"
}
### VIASH END

input_path = f"{meta['resources_dir']}/pancreas/dataset.h5ad"
input = sc.read_h5ad(input_path)

print(">> Running script as test for even")
output_path = "output.h5ad"
out = subprocess.check_output([
    meta["executable"],
    "--input", input_path,
    "--output", output_path,
    "--even",
    "--seed", "123",
    "--n_obs", "100",
    "--n_vars", "120"
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_path)

print(">> Check that test output fits expected API")
output = sc.read_h5ad(output_path)

assert output.n_obs <= 100
assert output.n_vars <= 120



print(">> Runing script as test for specific batch and celltype categories")
output2_path = "output.h5ad"

keep_features = list(input.var_names[:10])
out = subprocess.check_output([
    meta["executable"],
    "--input", input_path,
    "--keep_celltype_categories", "acinar:beta",
    "--keep_batch_categories", "celseq:inDrop4:smarter",
    "--keep_features", ":".join(keep_features),
    "--output", output_path,
    "--seed", "123"
]).decode("utf-8")
print(">> Checking whether file exists")
assert path.exists(output2_path)

print(">> Check that test output fits expected API")
output2 = sc.read_h5ad(output2_path)

assert output2.n_obs <= 500
assert output2.n_vars <= 500
assert set(keep_features).issubset(output2.var_names)
