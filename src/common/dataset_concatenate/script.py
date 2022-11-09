import scanpy as sc

###VIASH START
par = {
    "inputs": ["../resources/pancreas/toy_data.h5ad", "../resources/pancreas/toy_data.h5ad"],
    "output": "output.h5ad"
}
###VIASH END

adata_list = []
for i in par["inputs"]:
    print("Loading {}".format(i))
    adata_list.append(
        sc.read_h5ad(i)
    )

print("Concatenate anndatas")
adata = adata_list[0].concatenate(adata_list[1:])

print("Writing result file")
adata.write_h5ad(par["output"], compression="gzip")
