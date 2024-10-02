import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import anndata as ad
import scipy
from GPcounts.RNA_seq_GP import rna_seq_gp
import warnings
warnings.filterwarnings('ignore')


# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad',
    'n_features': 120
}
meta = {
    'functionality_name': 'GPcounts'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run GPcounts')

# Subset if required
if par['n_features']:
    adata = adata[:, :par['n_features']]

counts = adata.layers["counts"]
if scipy.sparse.issparse(counts): 
    counts = counts.todense()

# spatialx = [str(i) for i in adata.obsm['spatial'][:, 0]]
# spatialy = [str(i) for i in adata.obsm['spatial'][:, 1]]

# index_names = [i+'x'+j for i, j in zip(spatialx, spatialy)]
# Y = pd.DataFrame(data=counts, index=index_names, columns=adata.var.index)

# spatial_locations = pd.DataFrame(index=Y.index)
# spatial_locations['x'] = Y.index.str.split('x').str.get(0).map(float)
# spatial_locations['y'] = Y.index.str.split('x').str.get(1).map(float)

# spatial_locations['total_counts'] = Y.sum(1)

Y = pd.DataFrame(data=counts, 
                index=adata.obs_names, 
                columns=adata.var_names)

spatial_locations = pd.DataFrame(data=adata.obsm['spatial'], 
                                index=adata.obs_names, 
                                columns=['x', 'y'])
spatial_locations['total_counts'] = Y.sum(1)

Y = Y.loc[spatial_locations.index]
X = spatial_locations[['x', 'y']]

scales = []
for i in range(0, len(Y.columns)):
    model = smf.glm(formula="Y.iloc[:,i]~0+spatial_locations['total_counts']", data=Y,
                    family=sm.families.NegativeBinomial(sm.families.links.identity())).fit()
    res = model.params[0]*spatial_locations['total_counts']
    scales.append(res)
scalesdf = pd.DataFrame(scales)
scalesdf = scalesdf.T

Y = Y.T
X = X[['x', 'y']]

sparse = True
nb_scaled = True  # set the nb_scaled argument to True to pass the scale factors
gene_name = Y.index
likelihood = 'Negative_binomial'
gp_counts = rna_seq_gp(
    X, Y.loc[gene_name], sparse=sparse, M=250, scale=scalesdf, safe_mode=False)

log_likelihood_ratio = gp_counts.One_sample_test(likelihood)

df = gp_counts.calculate_FDR(log_likelihood_ratio)

# save results
df = df.loc[adata.var_names][['log_likelihood_ratio']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
