import os
import subprocess
import anndata as ad

input_expression ="https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5"
input_spatial = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_spatial.tar.gz"
dataset_id = "mouse_brain_coronal"
dataset_name = "Mouse Brain Coronal Section 1 (FFPE)"
dataset_url = "https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard"
dataset_summary = "Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set"
dataset_description = "CytAssist_FFPE_Mouse_Brain_Rep1 - Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set"
dataset_organism = "Mus musculus"
dataset = "dataset.h5ad"

print(">> Running script", flush=True)
out = subprocess.run(
    [
        meta['executable'],
        "--input_expression",  input_expression, 
        "--input_spatial", input_spatial, 
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
