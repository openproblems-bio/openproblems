import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from somde import SomNode
import scipy

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'SOMDE'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run SOMDE', flush=True)
counts = adata.layers["counts"]
if scipy.sparse.issparse(counts): 
    counts = counts.todense()
    
data = pd.DataFrame(
    counts, 
    columns=adata.var_names, 
    index=adata.obs_names
)

X = pd.DataFrame(adata.obsm["spatial"], 
                     index=adata.obs_names, 
                     columns=["x", "y"]).values.astype(np.float32)
    
som = SomNode(X, k=10)
ndf, ninfo = som.mtx(data.transpose())
nres = som.norm() 

df, SVnum = som.run()

# save results
df.set_index("g", inplace=True)
df = df.loc[adata.var_names][['FSV']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
