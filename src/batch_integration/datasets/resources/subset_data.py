import scanpy as sc

adata_file = 'src/common/data_loader/resources/pancreas.h5ad'
adata_out_file = 'src/batch_integration/datasets/resources/data_loader_pancreas.h5ad'
batch = 'tech'
label = 'celltype'

adata = sc.read(adata_file)
print(adata)

# observation subset
head_batches = adata.obs[batch].unique().tolist()[0:3]
head_labels = adata.obs[label].unique().tolist()[0:2]

# feature subset
g2m_file = 'src/batch_integration/datasets/resources/g2m_genes_tirosh_hm.txt'
s_file = 'src/batch_integration/datasets/resources/s_genes_tirosh_hm.txt'
g2m_genes = [x.strip() for x in open(g2m_file).readlines()]
s_genes = [x.strip() for x in open(s_file).readlines()]

all_genes = adata.var.index.tolist()
cc_genes = [x for x in g2m_genes + s_genes if x in all_genes]
head_genes = list(set(cc_genes + all_genes[:200]))

# subset adata
adata = adata[adata.obs[batch].isin(head_batches)]
adata = adata[adata.obs[label].isin(head_labels)]
adata = adata[:, head_genes]

print(adata)

adata.write(adata_out_file, compression='gzip')
