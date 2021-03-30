import scanpy as sc

def preprocess_scanpy(adata):

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')

    return adata
