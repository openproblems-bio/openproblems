import sys
import anndata as ad

## VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    "output_train": "train.h5ad",
    "output_test": "test.h5ad",
}
meta = {
    "functionality_name": "split_data",
    "config": "src/dimensionality_reduction/process_dataset/.config.vsh.yaml"
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
output_train = subset_anndata(adata, slot_info["output_train"])

print(">> Creating test data", flush=True)
output_test = subset_anndata(adata, slot_info["output_test"])

print(">> Writing", flush=True)
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
