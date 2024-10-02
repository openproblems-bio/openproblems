import anndata as ad
import numpy as np

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'input_solution': 'resources_test/spatially_variable_genes/mouse_brain_coronal/solution.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'random_ranking'
}
# VIASH END

print('Generate predictions', flush=True)
input_data = ad.read_h5ad(par['input_data'])

df = input_data.var[["feature_id"]]

np.random.seed(0)
df['pred_spatial_var_score'] = np.random.rand(len(df['feature_id']))

output = ad.AnnData(var=df,
                    uns={'dataset_id': input_data.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
