print("Loading dependencies", flush=True)
import scanpy as sc
import numpy as np

print("Reading input h5ad file", flush=True)
adata = sc.read_h5ad(par["input"])

print("Check parameters", flush=True)
new_shape = (adata.X.shape[0], 10)
adata.obsm["aligned"] = np.random.normal(0, 0.1, new_shape)
adata.obsm["mode2_aligned"] = np.random.normal(0, 0.1, new_shape)

print("Write output to file", flush=True)
adata.uns["method_id"] = "sample_method"
adata.write_h5ad(par["output"], compression = "gzip")
