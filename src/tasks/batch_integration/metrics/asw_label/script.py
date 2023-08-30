import anndata as ad
from scib.metrics import silhouette

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
input_solution.obsm["X_emb"] = input_integrated.obsm["X_emb"]

print('compute score', flush=True)
score = silhouette(
    input_solution,
    label_key='label',
    embed='X_emb'
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_solution.uns['dataset_id'],
        'normalization_id': input_solution.uns['normalization_id'],
        "method_id": input_integrated.uns['method_id'],
        "metric_ids": [meta['functionality_name']],
        "metric_values": [score]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
