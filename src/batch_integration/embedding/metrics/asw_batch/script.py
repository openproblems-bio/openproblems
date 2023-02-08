import pprint
import anndata as ad
from scib.metrics import silhouette_batch
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
    'config': 'src/batch_integration/embedding/metrics/asw_batch/config.vsh.yaml'
}
## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["output_type"]
integrated_embedding = config["functionality"]["info"]["integrated_embedding"]

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

print('compute score')
score = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='label',
    embed=integrated_embedding,
)

print("Create output AnnData object")
output = ad.AnnData(
    uns={
        "dataset_id": adata.uns['dataset_id'],
        "method_id": adata.uns['method_id'],
        "metric_ids": [ meta['functionality_name'] ],
        "metric_values": [ score ],
        "hvg": adata.uns['hvg'],
        "scaled": adata.uns['scaled'],
        "output_type": output_type,
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
