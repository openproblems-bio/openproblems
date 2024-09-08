import sys
import anndata as ad
import scanpy as sc
from scib.metrics.clustering import cluster_optimal_resolution
from scib.metrics import ari, nmi

## VIASH START
par = {
    'adata_integrated': 'resources_test/batch_integration/pancreas/integrated_graph.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsp='obsp', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('Run optimal Leiden clustering', flush=True)
cluster_optimal_resolution(
    adata=adata,
    label_key='label',
    cluster_key='cluster',
    cluster_function=sc.tl.leiden,
)

print('Compute ARI score', flush=True)
ari_score = ari(adata, cluster_key='cluster', label_key='label')

print('Compute NMI score', flush=True)
nmi_score = nmi(adata, cluster_key='cluster', label_key='label')

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        "method_id": adata.uns['method_id'],
        "metric_ids": [ "ari", "nmi" ],
        "metric_values": [ ari_score, nmi_score ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")