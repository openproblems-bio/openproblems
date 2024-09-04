import anndata as ad
import umap

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "reduced.h5ad",
    "n_comps": 2,
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

def diffusion_map(graph, n_comps, t, n_retries):
    import numpy as np
    import scipy.sparse.linalg

    diag_data = np.asarray(graph.sum(axis=0))
    identity = scipy.sparse.identity(graph.shape[0], dtype=np.float64)
    diag = scipy.sparse.spdiags(
        1.0 / np.sqrt(diag_data), 0, graph.shape[0], graph.shape[0]
    )
    laplacian = identity - diag * graph * diag
    num_lanczos_vectors = max(2 * n_comps + 1, int(np.sqrt(graph.shape[0])))
    try:
        eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
            laplacian,
            n_comps,
            which="SM",
            ncv=num_lanczos_vectors,
            tol=1e-4,
            v0=np.ones(laplacian.shape[0]),
            maxiter=graph.shape[0] * 5,
        )
        return (eigenvalues**t) * eigenvectors
    except scipy.sparse.linalg.ArpackNoConvergence:
        if n_retries > 0:
            # add some noise and try again
            graph_rand = graph.copy().tocoo()
            graph_rand.row = np.random.choice(
                graph_rand.shape[0], len(graph_rand.row), replace=True
            )
            graph_rand.data *= 0.01
            return diffusion_map(
                graph + graph_rand, n_comps, t, n_retries=n_retries - 1
            )
        else:
            raise

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

print("Create high dimensionally embedding with all features", flush=True)

n_comps = min(par["n_comps"], min(input.shape) - 2)

graph = umap.UMAP(transform_mode="graph").fit_transform(input.layers["normalized"])

X_emb = diffusion_map(graph, n_comps, t=par["t"], n_retries=par["n_retries"])


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