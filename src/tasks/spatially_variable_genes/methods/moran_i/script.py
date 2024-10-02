import warnings
warnings.filterwarnings('ignore')

import anndata as ad
import squidpy as sq

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad',
    'coord_type_moran_i': 'generic'
    
}
meta = {
    'functionality_name': 'moranI'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run moranI', flush=True)
sq.gr.spatial_neighbors(adata,
                        coord_type=par['coord_type_moran_i'],
                        delaunay=True)

sq.gr.spatial_autocorr(adata,
                       mode="moran",
                       layer='normalized',
                       n_perms=100,
                       genes=adata.var_names)

# save results
df = adata.uns["moranI"]
df = df.loc[adata.var_names][['I']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
