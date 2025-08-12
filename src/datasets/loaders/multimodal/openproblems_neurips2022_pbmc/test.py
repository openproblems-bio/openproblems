from os import path
import subprocess
import anndata as ad

# TODO: update once data is public

input_mod1 = "cite_rna_merged.h5ad" #change data set path after loading manually?
input_mod2 = "cite_prot_merged.h5ad" #change data set path after loading manually?
mod1 = "GEX"
mod2 = "ADT"

output_mod1_file = "output_mod1.h5ad"
output_mod2_file = "output_mod2.h5ad"

input_url_mod1 = "s3://openproblems-nextflow/datasets_private/neurips2022/cite_rna_merged.h5ad"
input_url_mod2 = "s3://openproblems-nextflow/datasets_private/neurips2022/cite_prot_merged.h5ad"

# download input
# print(">> Downloading input modality 1", flush=True)
# out = subprocess.run(
#   [
#     "aws s3 cp",
#     "-O", input_mod1,
#     input_url_mod1,
#   ],
#   stderr=subprocess.STDOUT
# )

# print(">> Downloading input modality 2", flush=True)
# out = subprocess.run(
#   [
#     "aws s3 cp",
#     "-O", input_mod2,
#     input_url_mod2,
#   ],
#   stderr=subprocess.STDOUT
# )


print(">> Running script", flush=True)
out = subprocess.run(
  [
    meta["executable"],
    "--input_mod1", input_mod1,
    "--input_mod2", input_mod2,
    "--mod1", mod1,
    "--mod2", mod2,
    "--output_mod1", output_mod1_file,
    "--output_mod2", output_mod2_file,
    "--dataset_id", "openproblems/neurips2021_bmmc",
    "--dataset_name", "Kaggle22 PBMC (CITE-seq)",
    "--dataset_url", "https://www.kaggle.com/competitions/open-problems-multimodal/data",
    "--dataset_reference", "Neurips22",
    "--dataset_summary", "Neurips22 competition dataset",
    "--dataset_description", "The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors.",
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
assert output_mod1.uns["dataset_name"] == "Kaggle22 PBMC (CITE-seq)", "Expected: Kaggle22 PBMC (CITE-seq) as value for dataset_name in mod 1 output uns"
assert output_mod1.uns["dataset_url"] == "https://www.kaggle.com/competitions/open-problems-multimodal/data", "Expected: https://www.kaggle.com/competitions/open-problems-multimodal/data as value for dataset_url in mod 1 output uns"
assert output_mod1.uns["dataset_reference"] == "Neurips22", "Expected: Neurips22 as value for dataset_reference in mod 1 output uns"
assert output_mod1.uns["dataset_summary"] == "Neurips22 competition dataset", "Expected: Neurips22 competition dataset as value for dataset_summary in mod 1 output uns"
assert output_mod1.uns["dataset_description"] == "The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors.", "Expected: The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors. as value for dataset_description in mod 1 output uns"


print(">> Check that output mod2 fits expected API", flush=True)
assert output_mod2.X is None, ".X is not None/empty in mod 2 output"
assert "counts" in output_mod2.layers, "'counts' not found in mod 2 output layers"
assert "cell_type" in output_mod2.obs.columns, "cell_type column not found in mod 2 output obs"
assert "batch" in output_mod2.obs.columns, "batch column not found in mod 2 output obs"
assert output_mod2.uns["dataset_name"] == "Kaggle22 PBMC (CITE-seq)", "Expected: Kaggle22 PBMC (CITE-seq) as value for dataset_name in mod 2 output uns"
assert output_mod2.uns["dataset_url"] == "https://www.kaggle.com/competitions/open-problems-multimodal/data", "Expected: https://www.kaggle.com/competitions/open-problems-multimodal/data as value for dataset_url in mod 2 output uns"
assert output_mod2.uns["dataset_reference"] == "Neurips22", "Expected: Neurips22 as value for dataset_reference in mod 2 output uns"
assert output_mod2.uns["dataset_summary"] == "Neurips22 competition dataset", "Expected: Neurips22 competition dataset as value for dataset_summary in mod 2 output uns"
assert output_mod2.uns["dataset_description"] == "The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors.", "Expected: The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors. as value for dataset_description in mod 2 output uns"