import h5py
import numpy as np
import anndata as ad
import subprocess
from os import path

input_train_path = meta["resources_dir"] + "/pancreas/train.h5ad"
output_path = "output.h5ad"

cmd = [
    meta['executable'],
    "--input_train", input_train_path,
    "--output", output_path
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

# helper functions interpreted from https://stackoverflow.com/questions/44883175/how-to-list-all-datasets-in-h5py-file
def descend_obj(obj,sep='\t'):
    """
    Iterate through groups in a HDF5 file and prints the groups and datasets names and datasets attributes
    """
    if type(obj) in [h5py._hl.group.Group,h5py._hl.files.File]:
        for key in obj.keys():
            print(f"{sep}{key}: {obj[key]}", flush=True)
            descend_obj(obj[key],sep=sep+'\t')
    elif type(obj)==h5py._hl.dataset.Dataset:
        for key in obj.attrs.keys():
            print(f"{sep}\t{key}: {obj.attrs[key]}", flush=True)

print(">> Reading h5ad files")
with h5py.File(input_train_path, 'r') as input_train:
    with h5py.File(output_path, 'r') as output:
        print("Input:", flush=True)
        descend_obj(input_train)
        print("Output:", flush=True)
        descend_obj(output)

        print(">> Checking whether predictions were added", flush=True)
        assert "denoised" in output["layers"].keys()
        assert meta['functionality_name'] == np.string_(output["uns/method_id"]).decode('utf-8')

        print("Checking whether data from input was copied properly to output")
        assert input_train.get("layers/counts").shape == output.get("layers/counts").shape
        assert input_train["uns/dataset_id"] == output["uns/dataset_id"]

print("All checks succeeded!")
