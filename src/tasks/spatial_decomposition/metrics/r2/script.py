import anndata as ad
import sklearn.metrics

## VIASH START
par = {
  'input_method': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/output.h5ad',
  'input_solution': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/solution.h5ad',
  'output': 'score.h5ad'
}
meta = {
  'functionality_name': 'r2'
}
## VIASH END

print('Reading input files', flush=True)
input_method = ad.read_h5ad(par['input_method'])
input_solution = ad.read_h5ad(par['input_solution'])

print('Compute metrics', flush=True)
prop_true = input_solution.obsm["proportions_true"]
prop_pred = input_method.obsm["proportions_pred"]
r2_score = sklearn.metrics.r2_score(
  prop_true, prop_pred, sample_weight=None, multioutput="uniform_average"
)

uns_metric_ids = [ 'r2' ]
uns_metric_values = [ r2_score ]

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

