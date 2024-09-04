import os
import subprocess
import anndata as ad

input_data ="https://zenodo.org/records/12784832/files/Slide-seqV2_stickels2020highly_stickels2021highly_SlideSeqV2_Mouse_Olfactory_bulb_Puck_200127_15_data_whole.h5ad?download=1"
dataset_id = "zenodo_spatial/mouse_olfactory_bulb_puck"
dataset_name = "mouse_olfactory_bulb_puck"
dataset_url = "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
dataset_summary = "Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2"
dataset_description = "Gene expression library of mouse olfactory bulk puck profiled using Slide-seq V2"
dataset_organism = "Mus musculus"
dataset = "dataset.h5ad"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta['executable'],
        "--input_data",  input_data, 
        "--dataset_id", dataset_id, 
        "--dataset_name", dataset_name, 
        "--dataset_url", dataset_url, 
        "--dataset_summary", dataset_summary, 
        "--dataset_description", dataset_description, 
        "--dataset_organism", dataset_organism, 
        "--dataset", dataset
    ],
    stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout, flush=True)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.", flush=True)
    exit(out.returncode)

print(">> Checking whether output file exists", flush=True)
assert os.path.exists(dataset), "Output does not exist"

print(">> Read output anndata", flush=True)
adata = ad.read_h5ad(dataset)

print(adata)

print(">> Check that output fits expected API", flush=True)
assert adata.X is None, "adata.X should be None/empty"
assert "counts" in adata.layers, "Counts layer not found in .layers"
assert adata.uns["dataset_id"] == dataset_id, f"Expected {dataset_id} as value"
assert adata.uns["dataset_name"] == dataset_name, f"Expected {dataset_name} as value"
assert adata.uns["dataset_url"] == dataset_url, f"Expected {dataset_url} as value"
assert adata.uns["dataset_summary"] == dataset_summary, f"Expected {dataset_summary} as value"
assert adata.uns["dataset_organism"] == dataset_organism, f"Expected {dataset_organism} as value"
assert 'spatial' in adata.obsm, "Spatial spot coordinates not found in .obsm"

print(">> All tests passed successfully", flush=True)
