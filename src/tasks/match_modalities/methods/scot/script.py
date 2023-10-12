import anndata as ad
import sys
sys.path.append("/opt/SCOT/src/")
import scotv1
import pandas as pd

# importing helper functions from common preprocessing.py file in resources dir
import sys


## VIASH START
par = {
  "input_mod1" : "resources_test/common/scicar_cell_lines/dataset_mod1.h5ad",
  "input_mod2" : "resources_test/common/scicar_cell_lines/dataset_mod2.h5ad",
  "output_mod1" : "integrated_mod1.h5ad",
  "output_mod2" : "integrated_mod2.h5ad",
  "balanced":False,
}
## VIASH END


print("Reading input h5ad file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])


print("Initialize SCOT", flush=True)
scot = scotv1.SCOT(adata_mod1.obsm["X_svd"], adata_mod2.obsm["X_svd"])

print("Call the unbalanced alignment", flush=True)
# From https://github.com/rsinghlab/SCOT/blob/master/examples/unbalanced_GW_SNAREseq.ipynb # noqa: 501
X_new_unbal, y_new_unbal = scot.align(
    k=50, e=1e-3, normalize=True
)


print("store output", flush=True)
adata_mod1.obsm["integrated"] = X_new_unbal
adata_mod2.obsm["integrated"] = y_new_unbal

print("Write output to file", flush=True)
adata_mod1.uns["method_id"] = meta["functionality_name"]
adata_mod2.uns["method_id"] = meta["functionality_name"]
adata_mod1.write_h5ad(par["output_mod1"], compression = "gzip")
adata_mod2.write_h5ad(par["output_mod2"], compression = "gzip")
