import anndata as ad
import squidpy as sq

### VIASH START
par = {
    "input": "resources_test/common/mouse_brain_coronal_section1/dataset.h5ad",
    "input_layer": "normalized",
    "output": "reference_dataset.h5ad",
    "num_features": 50,
    "coord_type_proc": "grid"
}
### VIASH END

print(">> Load data", flush=True)
adata = ad.read_h5ad(par['input'])

print(">> Run Moran's I spatial autocorrelation", flush=True)
sq.gr.spatial_neighbors(adata, 
                        coord_type=par['coord_type_proc'], 
                        delaunay=False)
sq.gr.spatial_autocorr(adata, 
                       layer="normalized",
                       mode="moran", 
                       n_perms=100, n_jobs=10, 
                       genes=adata.var_names)

n_svgs = par['num_features']
sel_genes = (
    adata.uns["moranI"]["I"].sort_values(ascending=False).head(n_svgs).index.tolist()
)

adata = adata[:, sel_genes]

print(">> Writing data", flush=True)
adata.write_h5ad(par['output'])

