import os
from os import path
import subprocess

import scanpy as sc
import pandas
import numpy as np

print(">> Running sample_dataset")
out = subprocess.check_output([
    "./sample_dataset",
    "--output", "output.h5ad"
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists("output.h5ad")

print(">> Check that output fits expected API")
adata = sc.read_h5ad("output.h5ad")
assert "mode2" in adata.obsm
assert "mode2_obs" in adata.uns
assert "mode2_var" in adata.uns
assert np.all(adata.obs.index == adata.uns["mode2_obs"])
assert len(adata.uns["mode2_var"]) == adata.obsm["mode2"].shape[1]

# check dataset id
assert "dataset_id" in adata.uns




print(">> Running sample_dataset with different args")
out = subprocess.check_output([
    "./sample_dataset",
    "--output", "output.h5ad",
    "--n_cells", "100",
    "--n_genes", "200"
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists("output.h5ad")

print(">> Check that output fits expected API")
adata = sc.read_h5ad("output.h5ad")
assert "mode2" in adata.obsm
assert "mode2_obs" in adata.uns
assert "mode2_var" in adata.uns
assert np.all(adata.obs.index == adata.uns["mode2_obs"])
assert len(adata.uns["mode2_var"]) == adata.obsm["mode2"].shape[1]

# check shape based on args
assert adata.shape == (100, 200)

# check dataset id
assert "dataset_id" in adata.uns



print(">> All tests passed successfully")
