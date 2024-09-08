import sys
import anndata as ad
from scib.metrics import hvg_overlap

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
adata_solution = read_anndata(
    par['input_solution'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)
adata_integrated = read_anndata(
    par['input_integrated'],
    X='layers/corrected_counts',
    obs='obs',
    var='var',
    uns='uns'
)

print('compute score', flush=True)
score = hvg_overlap(
    adata_solution,
    adata_integrated,
    batch_key="batch"
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata_solution.uns['dataset_id'],
        'normalization_id': adata_solution.uns['normalization_id'],
        "method_id": adata_integrated.uns['method_id'],
        "metric_ids": [meta['functionality_name']],
        "metric_values": [score]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
