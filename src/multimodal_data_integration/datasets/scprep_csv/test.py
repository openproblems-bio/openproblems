import os
from os import path
import subprocess

import scanpy as sc
import pandas
import numpy as np

import urllib.request

print(">> Downloading input file")
# need to download file manually for now; viash docker platform tries to auto-mount them
urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FCD8%5Fmerged%2DADT%5Fumi%2Ecsv%2Egz", "adt_umi.csv.gz")

print(">> Running scprep_csv")

out = subprocess.check_output([
    "./scprep_csv",
    "--id", "footest",
    "--input1", "adt_umi.csv.gz",
    "--input2", "adt_umi.csv.gz",
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

# since same file was used for both datasets
assert adata.shape == adata.obsm["mode2"].shape

# check dataset id
assert "dataset_id" in adata.uns
assert adata.uns["dataset_id"] == "footest"

print(">> All tests passed successfully")
