import sys
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'ouput': 'output.h5ad'
}

meta = { 
    'functionality': 'foo',
    'config': 'bar'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    obs='obs',
    obsm='obsm',
    uns='uns'
)


print('Run kNN...', flush=True)
sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.write_h5ad(par['output'], compression='gzip')