import scanpy as sc
import yaml

## VIASH START

par = {
    'input': 'resources_test/batch_integration/combat.h5ad',
    'ouput': 'output.h5ad'
}

meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}

## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["output_type"]

print('Read input', flush=True)
adata= sc.read_h5ad(par['input'])


print('Run PCA', flush=True)
adata.obsm['X_emb'] = sc.pp.pca(
    adata.layers["corrected_counts"],
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=False
)

print('Store outputs', flush=True)
adata.uns['output_type'] = output_type
adata.write_h5ad(par['output'], compression='gzip')