from os import path
import subprocess
import anndata as ad

input = "neurips2021_bmmc_cite.h5ad"
mod1 = "GEX"
mod2 = "ADT"

output_mod1_file = "output_mod1.h5ad"
output_mod2_file = "output_mod2.h5ad"

input_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz"

# download input
print(">> Downloading input", flush=True)
out = subprocess.run(
  [
    "wget",
    "-O", input + ".gz",
    input_url,
  ],
  stderr=subprocess.STDOUT
)
# unzip input
print(">> Unzipping input", flush=True)
out = subprocess.run(
  [
    "gunzip",
    input + ".gz",
  ],
  stderr=subprocess.STDOUT
)

print(">> Running script", flush=True)
out = subprocess.run(
  [
    meta["executable"],
    "--input", input,
    "--mod1", mod1,
    "--mod2", mod2,
    "--output_mod1", output_mod1_file,
    "--output_mod2", output_mod2_file,
    "--dataset_id", "openproblems/neurips2021_bmmc",
    "--dataset_name", "BMMC (Multiome)",
    "--dataset_url", "http://foo.org",
    "--dataset_reference", "foo2000bar",
    "--dataset_summary", "A short summary.",
    "--dataset_description", "A couple of paragraphs worth of text.",
    "--dataset_organism", "homo_sapiens",
  ],
  stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout, flush=True)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.", flush=True)
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
assert "cell_type" in output_mod1.obs.columns, "cell_type column not found in mod 1 output obs"
assert "batch" in output_mod1.obs.columns, "batch column not found in mod 1 output obs"
assert output_mod1.uns["dataset_name"] == "BMMC (Multiome)", "Expected: Pancreas as value for dataset_name in mod 1 output uns"
assert output_mod1.uns["dataset_url"] == "http://foo.org", "Expected: http://foo.org as value for dataset_url in mod 1 output uns"
assert output_mod1.uns["dataset_reference"] == "foo2000bar", "Expected: foo2000bar as value for dataset_reference in mod 1 output uns"
assert output_mod1.uns["dataset_summary"] == "A short summary.", "Expected: A short summary. as value for dataset_summary in mod 1 output uns"
assert output_mod1.uns["dataset_description"] == "A couple of paragraphs worth of text.", "Expected: A couple of paragraphs worth of text. as value for dataset_description in mod 1 output uns"


print(">> Check that output mod2 fits expected API", flush=True)
assert output_mod2.X is None, ".X is not None/empty in mod 2 output"
assert "counts" in output_mod2.layers, "'counts' not found in mod 2 output layers"
assert "cell_type" in output_mod2.obs.columns, "cell_type column not found in mod 2 output obs"
assert "batch" in output_mod2.obs.columns, "batch column not found in mod 2 output obs"
assert output_mod2.uns["dataset_name"] == "BMMC (Multiome)", "Expected: Pancreas as value for dataset_name in mod 2 output uns"
assert output_mod2.uns["dataset_url"] == "http://foo.org", "Expected: http://foo.org as value for dataset_url in mod 2 output uns"
assert output_mod2.uns["dataset_reference"] == "foo2000bar", "Expected: foo2000bar as value for dataset_reference in mod 2 output uns"
assert output_mod2.uns["dataset_summary"] == "A short summary.", "Expected: A short summary. as value for dataset_summary in mod 2 output uns"
assert output_mod2.uns["dataset_description"] == "A couple of paragraphs worth of text.", "Expected: A couple of paragraphs worth of text. as value for dataset_description in mod 2 output uns"