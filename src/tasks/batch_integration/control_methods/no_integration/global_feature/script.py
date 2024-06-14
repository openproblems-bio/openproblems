import scanpy as sc

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
}

meta = { 
    'functionality': 'foo',
    'config': 'bar',
    "resources_dir": "src/tasks/batch_integration/control_methods/"
}

## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

# no processing, subset matrix to highly variable genes
adata_hvg = adata[:, adata.var["hvg"]].copy()
adata.layers['corrected_counts'] = adata_hvg.layers["normalized"].copy()

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
