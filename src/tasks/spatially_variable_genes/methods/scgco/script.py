import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import anndata as ad
import numpy as np
import scipy
import sys
sys.path.append("/opt/scGCO")

from scGCO_simple import *

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal_section1/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'scGCO'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

counts = adata.layers["counts"]
if scipy.sparse.issparse(counts): 
    counts = counts.todense()

data = pd.DataFrame(
    counts,
    columns=adata.var_names,
    index=adata.obs_names
)

print('Run scGCO', flush=True)
data_norm = normalize_count_cellranger(data)

exp = data.iloc[:, 0]
locs = adata.obsm['spatial'].copy()

print('Create graph with weight', flush=True)
cellGraph = create_graph_with_weight(locs, exp)
gmmDict = gmm_model(data_norm)

print('Identify spatial genes', flush=True)
df = identify_spatial_genes(locs, data_norm, cellGraph, gmmDict)

# save results
df = df.loc[adata.var_names][['fdr']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

# Transform the values via -log10 to make sure a bigger score represents a 
# higher spatial variation
df['pred_spatial_var_score'] = -np.log10(df['pred_spatial_var_score'].tolist())

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
