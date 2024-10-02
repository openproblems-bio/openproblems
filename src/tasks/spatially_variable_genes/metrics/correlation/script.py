import anndata as ad
import pandas as pd

## VIASH START
par = {
  'input_method': 'resources_test/spatially_variable_genes/mouse_brain_coronal/output.h5ad',
  'input_solution': 'resources_test/spatially_variable_genes/mouse_brain_coronal/solution.h5ad',
  'output': 'score.h5ad'
}
meta = {
  'functionality_name': 'correlation'
}
## VIASH END

print('Reading input files', flush=True)
input_method = ad.read_h5ad(par['input_method'])
input_solution = ad.read_h5ad(par['input_solution'])

print('Compute metrics', flush=True)
df = pd.merge(input_method.var, input_solution.var, how='left', on='feature_id')
groupby = df.groupby('orig_feature_name', observed=True)
corr = groupby.apply(lambda x: x['pred_spatial_var_score'].corr(x['true_spatial_var_score'], method='kendall'))

uns_metric_ids = [ 'correlation' ]
uns_metric_values = [ corr.mean() ]

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'dataset_id': input_method.uns['dataset_id'],
    'method_id': input_method.uns['method_id'],
    'metric_ids': uns_metric_ids,
    'metric_values': uns_metric_values
  }
)
output.write_h5ad(par['output'], compression='gzip')

