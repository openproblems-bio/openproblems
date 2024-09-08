import sys
import anndata as ad
from scib.metrics import cell_cycle
import numpy as np

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END
sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata_solution = read_anndata(
    par['input_solution'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)
adata_integrated = read_anndata(
    par['input_integrated'],
    obs='obs',
    obsm='obsm',
    uns='uns'
)

print('Use gene symbols for features', flush=True)
adata_solution.var_names = adata_solution.var['feature_name']

translator = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
}

print('Compute score', flush=True)
if adata_solution.uns['dataset_organism'] not in translator:
    score = np.nan
else:
    organism = translator[adata_solution.uns['dataset_organism']]
    score = cell_cycle(
        adata_solution,
        adata_integrated,
        batch_key='batch',
        embed='X_emb',
        organism=organism,
    )

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata_solution.uns['dataset_id'],
        'normalization_id': adata_solution.uns['normalization_id'],
        'method_id': adata_integrated.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ]
    }
)


print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
