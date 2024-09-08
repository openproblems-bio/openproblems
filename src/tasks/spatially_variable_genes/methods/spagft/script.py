import anndata as ad
import SpaGFT as spg

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal_section1/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'SpaGFT'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run SpaGFT', flush=True)

adata.X = adata.layers['normalized'].copy()

adata.obs.loc[:, ['array_row', 'array_col']] = adata.obsm['spatial']

(ratio_low, ratio_high) = spg.gft.determine_frequency_ratio(adata, ratio_neighbors=1)

df = spg.detect_svg(adata,
                    spatial_info=['array_row', 'array_col'],
                    ratio_low_freq=ratio_low,
                    ratio_high_freq=ratio_high,
                    ratio_neighbors=1,
                    filter_peaks=True,
                    S=6)


# save results
df = df.loc[adata.var_names][['gft_score']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
