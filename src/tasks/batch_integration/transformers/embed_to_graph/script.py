import yaml
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'ouput': 'output.h5ad'
}
## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

print('Run kNN', flush=True)
sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.write_h5ad(par['output'], compression='gzip')