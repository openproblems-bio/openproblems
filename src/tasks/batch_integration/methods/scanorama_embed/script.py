import anndata as ad
import scanorama

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

# based on scib
# -> https://github.com/theislab/scib/blob/59ae6eee5e611d9d3db067685ec96c28804e9127/scib/utils.py#L51C1-L72C62
def merge_adata(*adata_list, **kwargs):
    """Merge adatas from list while remove duplicated ``obs`` and ``var`` columns

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """

    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(
                    f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`."
                )
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return ad.AnnData.concatenate(*adata_list, **kwargs)


print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par['n_hvg']:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:par['n_hvg']]
    adata = adata[:, idx].copy()

print('Run scanorama', flush=True)
adata.X = adata.layers['normalized']
split = []
batch_categories = adata.obs['batch'].cat.categories
for i in batch_categories:
    split.append(adata[adata.obs['batch'] == i].copy())
corrected = scanorama.correct_scanpy(split, return_dimred=True)
corrected = merge_adata(*corrected, batch_key='batch', batch_categories=batch_categories, index_unique=None)

print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    },
    obsm={
        'X_emb': corrected.obsm["X_scanorama"],
    }
)

print("Write output to file", flush=True)
output.write(par['output'], compression='gzip')
