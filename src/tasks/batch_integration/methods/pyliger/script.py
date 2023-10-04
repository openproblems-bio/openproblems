import anndata as ad
import numpy as np
import pyliger

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'pyliger'
}
## VIASH END

print('>> Read input', flush=True)
adata = ad.read_h5ad(par['input'])

print('>> Prepare data', flush=True)
adata_per_batch = []
for batch in adata.obs['batch'].unique():
  adb = adata[adata.obs['batch'] == batch].copy()

  # move counts
  adb.X = adb.layers['counts']
  del adb.layers['counts']

  # move normalized data
  adb.layers["norm_data"] = adb.layers["normalized"]
  del adb.layers["normalized"]
  
  # save row sum and sum of squares for further use
  norm_sum = np.ravel(np.sum(adb.layers["norm_data"], axis=0))
  norm_sum_sq = np.ravel(np.sum(adb.layers["norm_data"].power(2), axis=0))
  adb.var["norm_sum"] = norm_sum
  adb.var["norm_sum_sq"] = norm_sum_sq
  adb.var["norm_mean"] = norm_sum / adb.shape[0]

  # set more metadata
  adb.obs.index.name = 'cell_barcode'
  adb.var.index.name = 'gene_id'
  adb.uns['sample_name'] = batch

  # append to list
  adata_per_batch.append(adb)

print('Create liger object', flush=True)
lobj = pyliger.create_liger(
  adata_per_batch,
  remove_missing=False
)

# do not select genes
lobj.var_genes = adata.var_names

print('>> Scaling', flush=True)
pyliger.scale_not_center(lobj, remove_missing=False)

print('>> Optimize ALS', flush=True)
pyliger.optimize_ALS(lobj, k=20)

print('>> Quantile normalization', flush=True)
pyliger.quantile_norm(lobj)

print('>> Concatenate outputs', flush=True)
ad_out = ad.concat(lobj.adata_list)

print('>> Store output', flush=True)
adata.obsm['X_emb'] = ad_out[adata.obs_names, :].obsm['H_norm']
adata.uns['method_id'] = meta['functionality_name']

print("Write output to disk", flush=True)
adata.write_h5ad(par['output'], compression='gzip')
