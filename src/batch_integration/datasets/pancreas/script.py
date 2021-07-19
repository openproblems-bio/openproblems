## VIASH START
import os
print(os.getcwd())
par = {
    'adata': './src/batch_integration/resources/data_loader_pancreas.h5ad',
    'label': 'celltype',
    'batch': 'tech',
    'hvgs': 2000,
    'output': 'adata_out.h5ad',
    'debug': True
}
## VIASH END

print('Importing libraries')
import scanpy as sc
import scprep


def log_scran_pooling(adata):
    """Normalize data with scran via rpy2."""
    _scran = scprep.run.RFunction(
        setup="library('scran')",
        args="sce, min.mean=0.1",
        body="""
        sce <- computeSumFactors(
            sce, min.mean=min.mean,
            assay.type="X"
        )
        sizeFactors(sce)
        """,
    )
    adata.obs["size_factors"] = _scran(adata)
    adata.X = scprep.utils.matrix_vector_elementwise_multiply(
        adata.X, adata.obs["size_factors"], axis=0
    )
    sc.pp.log1p(adata)


if par['debug']:
    import pprint

    pprint.pprint(par)

adata_file = par['adata']
label = par['label']
batch = par['batch']
hvgs = par['hvgs']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)

# Rename columns
adata.obs['label'] = adata.obs[label]
adata.obs['batch'] = adata.obs[batch]
adata.layers['counts'] = adata.X

print(f'Select {hvgs} highly variable genes')
if adata.n_obs > hvgs:
    sc.pp.subsample(adata, n_obs=hvgs)

print('Normalisation with scran')
# only if "lognorm" exists in adata.layers
log_scran_pooling(adata)
adata.layers['logcounts'] = adata.X

print('Transformation: PCA')
sc.tl.pca(
    adata,
    svd_solver='arpack',
    return_info=True,
)
adata.obsm['X_uni'] = adata.obsm['X_pca']

print('Transformation: kNN')
sc.pp.neighbors(adata, use_rep='X_uni', key_added='uni')

print('Writing adata to file')
adata.write(output, compression='gzip')
