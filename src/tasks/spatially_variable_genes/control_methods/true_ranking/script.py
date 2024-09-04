import anndata as ad

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal_section1/dataset.h5ad',
    'input_solution': 'resources_test/spatially_variable_genes/mouse_brain_coronal_section1/solution.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'true_ranking'
}
# VIASH END

print('Generate predictions', flush=True)
input_solution = ad.read_h5ad(par['input_solution'])

df = input_solution.var[["feature_id", "true_spatial_var_score"]]
df.rename(columns={'true_spatial_var_score': 'pred_spatial_var_score'}, inplace=True)

output = ad.AnnData(var=df,
                    uns={'dataset_id': input_solution.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
