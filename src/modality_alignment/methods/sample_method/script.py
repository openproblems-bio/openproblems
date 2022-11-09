print("Loading dependencies")
import scanpy as sc
import numpy as np

print("Reading input h5ad file")
adata = sc.read_h5ad(par["input"])

print("Check parameters")
new_shape = (adata.X.shape[0], 10)
adata.obsm["aligned"] = np.random.normal(0, 0.1, new_shape)
adata.obsm["mode2_aligned"] = np.random.normal(0, 0.1, new_shape)

print("Write output to file")
adata.uns["method_id"] = "sample_method"
adata.write_h5ad(par["output"], compression = "gzip")
