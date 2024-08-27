import sys
import anndata as ad
from scib.metrics import kBET

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('compute score', flush=True)
score = kBET(
    adata,
    batch_key="batch",
    label_key="label",
    type_="embed",
    embed="X_emb",
    scaled=True,
    verbose=False,
)
print(score, flush=True)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')