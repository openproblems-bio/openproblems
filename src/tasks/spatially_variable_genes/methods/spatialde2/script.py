import scanpy as sc
import anndata as ad
import SpatialDE as sd
import NaiveDE
import warnings
warnings.filterwarnings("ignore")


# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'spatialDE2'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

# run SpatialDE2
print('Run spatialDE2', flush=True)
adata.X = adata.layers['counts'].copy()
sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[10])

counts = sc.get.obs_df(adata,
                       keys=list(adata.var_names),
                       use_raw=False,
                       layer='counts')

total_counts = sc.get.obs_df(adata, keys=["total_counts"])
norm_expr = NaiveDE.stabilize(counts.T).T
adata.X = NaiveDE.regress_out(
    total_counts, norm_expr.T, "np.log(total_counts)").T

# run SpatialDE2
df = sd.fit(adata, normalized=True, control=None)
df.set_index("gene", inplace=True)

# save results
df = df.loc[adata.var_names][['FSV']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
