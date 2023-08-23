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
adata = ad.read_h5ad(par['input_integrated'])

print('recompute kNN graph..', flush=True)
output_type = adata.uns['output_type']
adata_tmp = recompute_knn(
    adata,
    type_=output_type,
    use_rep= "X_emb" if output_type == 'embed' else "X_pca",
)

print('compute iLISI score...', flush=True)
ilisi_scores = lisi_graph_py(
    adata=adata_tmp,
    obs_key='batch',
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
ilisi = np.nanmedian(ilisi_scores)
ilisi = (ilisi - 1) / (adata.obs['batch'].nunique() - 1)

print('compute cLISI scores...', flush=True)
clisi_scores = lisi_graph_py(
    adata=adata_tmp,
    obs_key='label',
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
clisi = np.nanmedian(clisi_scores)
nlabs = adata.obs['label'].nunique()
clisi = (nlabs - clisi) / (nlabs - 1)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ 'ilisi_graph', 'clisi_graph' ],
        'metric_values': [ ilisi, clisi ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
