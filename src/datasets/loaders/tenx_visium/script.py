import subprocess
import squidpy as sq
import tempfile
import scanpy as sc

## VIASH START
par = {
  "input_expression": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5",
  "input_spatial": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_spatial.tar.gz",
  "dataset_id": "tenx_visium/mouse_brain_coronal_section1_visium",
  "dataset_name": "Mouse Brain Coronal Section 1 (FFPE)",
  "dataset_url": "https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard",
  "dataset_summary": "Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set",
  "dataset_organism": "Mus musculus",
  "dataset": "dataset.h5ad",
  "spot_filter_min_genes": 200,
  "gene_filter_min_spots": 50,
  "remove_mitochondrial": True
}
meta = {
  "functionality_name": "tenx_visium"
}
## VIASH END

print(f"Downloading data", flush=True)
with tempfile.TemporaryDirectory() as tempdir:
  input_exp = "feature_bc_matrix.h5"
  input_sp = "image_data.tar.gz"
  epx_data = subprocess.run(["wget", "-O", f"{tempdir}/{input_exp}", par['input_expression']], stderr=subprocess.STDOUT)
  sp_data = subprocess.run(["wget", "-O", f"{tempdir}/{input_sp}", par['input_spatial']], stderr=subprocess.STDOUT)
  extract_spatial = subprocess.run(["tar", "-xzf", f"{tempdir}/{input_sp}", "-C", tempdir], stderr=subprocess.STDOUT)

  # Read visium data and create anndata object
  adata = sq.read.visium(path=tempdir, counts_file=input_exp)

# Make variable names unique
adata.var_names_make_unique()

sc.pp.calculate_qc_metrics(adata, inplace=True)

print("Filtering spots or genes")
t0 = adata.shape
# remove cells with few counts
if par["spot_filter_min_counts"]:
  sc.pp.filter_cells(adata, min_counts=par["spot_filter_min_counts"], inplace=True)
# remove cells with few genes 
if par["spot_filter_min_genes"]:
  sc.pp.filter_cells(adata, min_genes=par["spot_filter_min_genes"], inplace=True)
# remove genes that have few counts
if par["gene_filter_min_counts"]:
  sc.pp.filter_genes(adata, min_counts=par["gene_filter_min_counts"], inplace=True)
# remove genes that are found in few cells
if par["gene_filter_min_spots"]:
  sc.pp.filter_genes(adata, min_cells=par["gene_filter_min_spots"], inplace=True)
t1 = adata.shape
print(f"Removed {t0[0] - t1[0]} cells and {(t0[1] - t1[1])} genes.")

if par["remove_mitochondrial"]:
  print("Removing mitochondrial genes")
  non_mito_genes_list = [name for name in adata.var_names if not (name.startswith('MT-') or name.startswith('mt-'))]
  adata = adata[:, non_mito_genes_list]
  
# Rename .var columns
adata.var['feature_name'] = adata.var_names
adata.var.set_index(adata.var['gene_ids'], inplace=True)
adata.var.rename(columns={"gene_ids": "feature_id"}, inplace=True)

# Move counts to .layers
print("Add metadata to uns", flush=True)
adata.layers["counts"] = adata.X
adata.X = None

# Add metadata
print("Add metadata to uns", flush=True)
metadata_fields = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for key in metadata_fields:
  if key in par:
    print(f"Setting .uns['{key}']", flush=True)
    adata.uns[key] = par[key]

print("Writing adata to file", flush=True)
adata.write_h5ad(par["dataset"], compression="gzip")