import anndata as ad

## VIASH START
par = {
    "inputs": ["resources_test/common/pancreas/dataset.h5ad", "resources_test/common/pancreas/dataset.h5ad"],
    "output": "output.h5ad"
}
## VIASH END

adata_list = []
for input in par["inputs"]:
    print("Loading {}".format(input))
    adata_list.append(ad.read_h5ad(input))

print("Concatenate anndatas")
adata = ad.concat(adata_list)

print("Concatenated anndata: ", adata)

print("Writing result file")
adata.write_h5ad(par["output"], compression="gzip")
