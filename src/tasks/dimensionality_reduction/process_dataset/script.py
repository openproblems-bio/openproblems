import sys
import anndata as ad

## VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    "output_dataset": "train.h5ad",
    "output_solution": "test.h5ad",
}
meta = {
    "functionality_name": "split_data",
    "config": "src/tasks/dimensionality_reduction/process_dataset/.config.vsh.yaml"
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_anndata import read_config_slots_info, subset_anndata

print(">> Load Data", flush=True)
adata = ad.read_h5ad(par["input"])

print(">> Figuring out which data needs to be copied to which output file", flush=True)
slot_info = read_config_slots_info(meta["config"])

print(">> Creating train data", flush=True)
output_dataset = subset_anndata(adata, slot_info["output_dataset"])

print(">> Creating test data", flush=True)
output_solution = subset_anndata(adata, slot_info["output_solution"])

print(">> Writing", flush=True)
output_dataset.write_h5ad(par["output_dataset"])
output_solution.write_h5ad(par["output_solution"])
