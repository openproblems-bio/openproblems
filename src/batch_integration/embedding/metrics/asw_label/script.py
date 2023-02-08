import pprint
import anndata as ad
from scib.metrics import silhouette
import yaml

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/processed.h5ad',
    'output': 'output.h5ad',
    'hvg': False,
    'scaling': False
}

meta = {
    'functionality_name': 'foo',
    'config': 'src/batch_integration/embedding/metrics/asw_label/config.vsh.yaml'
}
## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["output_type"]
integrated_embedding = config["functionality"]["info"]["integrated_embedding"]

adata_file = par['input']
output = par['output']

print('Read input', flush=True)
adata = ad.read_h5ad(adata_file)
name = adata.uns['dataset_id']

print('compute score')
score = silhouette(
    adata,
    group_key='label',
    embed=integrated_embedding
)

print("Create output AnnData object")
output = ad.AnnData(
    uns={
        "dataset_id": adata.uns['dataset_id'],
        "method_id": adata.uns['method_id'],
        "hvg": adata.uns['hvg'],
        "scaled": adata.uns['scaled'],
        "output_type": output_type,
        "metric_ids": meta['functionality_name'],
        "metric_value": score
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
