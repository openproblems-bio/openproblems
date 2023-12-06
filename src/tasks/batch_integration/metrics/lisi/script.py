import numpy as np
import anndata as ad
from scib.metrics.lisi import recompute_knn, lisi_graph_py

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print('Read input', flush=True)
input_solution = ad.read_h5ad(par['input_solution'])
input_integrated = ad.read_h5ad(par['input_integrated'])

input_solution.obsp["connectivities"] = input_integrated.obsp["connectivities"]
input_solution.obsp["distances"] = input_integrated.obsp["distances"]

# TODO: if we don't copy neighbors over, the metric doesn't work
input_solution.uns["neighbors"] = input_integrated.uns["neighbors"]

print('compute iLISI score...', flush=True)
ilisi_scores = lisi_graph_py(
    adata=input_solution,
    obs_key='batch',
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
ilisi = np.nanmedian(ilisi_scores)
ilisi = (ilisi - 1) / (input_solution.obs['batch'].nunique() - 1)

print('compute cLISI scores...', flush=True)
clisi_scores = lisi_graph_py(
    adata=input_solution,
    obs_key='label',
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
clisi = np.nanmedian(clisi_scores)
nlabs = input_solution.obs['label'].nunique()
clisi = (nlabs - clisi) / (nlabs - 1)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': input_solution.uns['dataset_id'],
        'normalization_id': input_solution.uns['normalization_id'],
        'method_id': input_integrated.uns['method_id'],
        'metric_ids': [ 'ilisi', 'clisi' ],
        'metric_values': [ ilisi, clisi ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
