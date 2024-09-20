import anndata as ad
import squidpy as sq

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad',
    'coord_type_sepal': 'grid',
    'max_neighs_sepal': 6,
}
meta = {
    'functionality_name': 'Sepal'
}
# VIASH END

print('Generate predictions', flush=True)
adata = ad.read_h5ad(par['input_data'])

sq.gr.spatial_neighbors(adata,
                        coord_type=par['coord_type_sepal'],
                        delaunay=False)

sq.gr.sepal(adata, 
            layer='normalized',
            max_neighs=par['max_neighs_sepal'], 
            genes=adata.var_names,
            n_jobs=1)

# save results
df = adata.uns["sepal_score"]
df = df.loc[adata.var_names][['sepal_score']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
