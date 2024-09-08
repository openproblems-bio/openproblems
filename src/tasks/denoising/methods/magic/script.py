import anndata as ad
import numpy as np
import scprep
from magic import MAGIC
import scipy


## VIASH START
par = {
    "input_train": "resources_test/denoising/pancreas/train.h5ad",
    "output": "output_magic.h5ad",
    "solver": "exact",
    "norm": "sqrt",
    "decay": 1,
    "t": 3,
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

print("Load data", flush=True)
input_train = ad.read_h5ad(par["input_train"], backed="r")

print("Set normalization method", flush=True)
if par["norm"] == "sqrt":
    norm_fn = np.sqrt
    denorm_fn = np.square
elif par["norm"] == "log":
    norm_fn = np.log1p
    denorm_fn = np.expm1
else:
    raise ValueError("Unknown normalization method: " + par["norm"] + ".")

print("Remove unneeded data", flush=True)
X = input_train.layers["counts"]

# Create output AnnData for later use
output = ad.AnnData(
    obs=input_train.obs[[]],
    var=input_train.var[[]],
    uns={
        "dataset_id": input_train.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    }
)

del input_train

print("Normalize data", flush=True)
X, libsize = scprep.normalize.library_size_normalize(
    X,
    rescale=1,
    return_library_size=True
)
X = scprep.utils.matrix_transform(X, norm_fn)

print("Run MAGIC", flush=True)
magic = MAGIC(
    solver=par["solver"],
    decay=par["decay"],
    t=par["t"],
    verbose=False,
)
X = magic.fit_transform(X, genes="all_genes")

print("Denormalizing data", flush=True)
X = scprep.utils.matrix_transform(X, denorm_fn)
X = scprep.utils.matrix_vector_elementwise_multiply(X, libsize, axis=0)

print("Create output AnnData", flush=True)
output.layers["denoised"] = X

print("Write Data", flush=True)
output.write_h5ad(par["output"], compression="gzip")

