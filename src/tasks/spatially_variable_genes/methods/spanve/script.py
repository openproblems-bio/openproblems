import anndata as ad
from Spanve import Spanve

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'Spanve'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run Spanve', flush=True)
adata.X = adata.layers['counts']
spanve = Spanve(adata)
spanve.fit(verbose=False)

# save results
df = spanve.result_df
df = df.loc[adata.var_names][['ent']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
