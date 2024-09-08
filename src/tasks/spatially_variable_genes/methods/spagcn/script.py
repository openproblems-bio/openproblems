import anndata as ad
import SpaGCN as spg
import pandas as pd
import numpy as np
import scanpy as sc
import random
import torch

# VIASH START
par = {
    'input_data': 'resources_test/spatially_variable_genes/mouse_brain_coronal_section1/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'SpaGCN'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

# run normalization
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata=adata)
sc.pp.log1p(adata)

print('Run SpaGCN', flush=True)
random_seed = 100

# Set seed
random.seed(random_seed)
torch.manual_seed(random_seed)
np.random.seed(random_seed)

p = 0.5
min_in_group_fraction = 0
min_in_out_group_ratio = 0
min_fold_change = 0


adj = spg.calculate_adj_matrix(
    x=adata.obsm["spatial"][:, 0],
    y=adata.obsm["spatial"][:, 1],
    histology=False
)
l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

clf = spg.SpaGCN()
clf.set_l(l)

# Run
clf.train(
    adata,
    adj,
    init_spa=True,
    init="louvain",
    res=0.5,
    tol=5e-3,
    lr=0.05,
    max_epochs=200,
)

y_pred, prob = clf.predict()
adata.obs["pred"] = y_pred
de_genes_all = list()
n_clusters = len(adata.obs["pred"].unique())

# identify DE genes
for target in range(n_clusters):
    print(f"target: {target}")
    start, end = np.quantile(adj[adj != 0], q=0.001), np.quantile(
        adj[adj != 0], q=0.1
    )
    r = spg.search_radius(
        target_cluster=target,
        cell_id=adata.obs.index.tolist(),
        x=adata.obsm["spatial"][:, 0],
        y=adata.obsm["spatial"][:, 1],
        pred=adata.obs["pred"].tolist(),
        start=start,
        end=end,
        num_min=10,
        num_max=14,
        max_run=100,
    )

    try:
        nbr_domians = spg.find_neighbor_clusters(
            target_cluster=target,
            cell_id=adata.obs.index.tolist(),
            x=adata.obsm["spatial"][:, 0],
            y=adata.obsm["spatial"][:, 1],
            pred=adata.obs["pred"].tolist(),
            radius=r,
            ratio=0,
        )

        de_genes_info = spg.rank_genes_groups(
            input_adata=adata,
            target_cluster=target,
            nbr_list=nbr_domians,
            label_col="pred",
            adj_nbr=True,
            log=True,
        )
        de_genes_all.append(de_genes_info)
    except (RuntimeError, TypeError, NameError):
        pass

if len(de_genes_all) == 0:
    df = adata.var
    df['pvals_adj'] = np.random.random(adata.n_vars)
else:
    df_res = pd.concat(de_genes_all)
    df_res = df_res.groupby(["genes"]).min()
    df_res = df_res.loc[adata.var_names]
    df = pd.concat([df_res, adata.var], axis=1)

# save results
df = df.loc[adata.var_names][['pvals_adj']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

# reverse it to make sure a bigger score represents a higher spatial variation
df['pred_spatial_var_score'] = -np.log10(df['pred_spatial_var_score'])

output = ad.AnnData(var=df,
                    uns={'dataset_id': adata.uns['dataset_id'],
                         'method_id': meta['functionality_name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
