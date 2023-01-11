import anndata as ad
import torch
from neuralee.embedding import NeuralEE
from neuralee.dataset import GeneExpressionDataset

# todo: allow gpu
device = torch.device("cpu")

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_hvg": 1000,
    "n_iter": 10,
    "normalize": True
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

if par["normalize"]:
    print("Performing own normalization", flush=True)
    # perform own normalization based on the "recommended" preprocessing taken from example notebooks, e.g.:
    # https://github.com/HiBearME/NeuralEE/blob/master/tests/notebooks/retina_dataset.ipynb
    dataset = GeneExpressionDataset(input.layers["counts"])
    dataset.log_shift()
    if par["n_hvg"]:
        dataset.subsample_genes(par["n_hvg"])
    dataset.standardscale()

else:
    X_mat = input.layers["normalized"]

    if par["n_hvg"]:
        print(f"Select top {par['n_hvg']} high variable genes", flush=True)
        idx = input.var["hvg_score"].to_numpy().argsort()[-par["n_hvg"]:]
        X_mat = X_mat[:, idx]
    
    print("Using pre-normalized data", flush=True)
    dataset = GeneExpressionDataset(X_mat)


# estimate the affinity matrix
batch_size = min(1000, input.n_obs)
print(f"Use {batch_size} cells as batch to estimate the affinity matrix", flush=True)
dataset.affinity_split(N_small=batch_size)

print("Create NeuralEE object", flush=True)
NEE = NeuralEE(dataset, d=2, device=device)
fine_tune_kwargs = dict(verbose=False)

if par["n_iter"]:
    fine_tune_kwargs["maxit"] = par["n_iter"]

print("Run NeuralEE", flush=True)
res = NEE.fine_tune(**fine_tune_kwargs)

X_emb = res["X"].detach().cpu().numpy()

print("Create output AnnData", flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    obsm={
        "X_emb": X_emb
    },
    uns={
        "dataset_id": input.uns["dataset_id"],
        "normalization_id": input.uns["normalization_id"],
        "method_id": meta["functionality_name"]
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")