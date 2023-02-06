from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.utils import check_version


@method(
    method_name="Diffusion maps",
    paper_reference="coifman2006diffusion",
    paper_name="Diffusion maps",
    paper_year=2006,
    code_url="https://github.com/openproblems-bio/openproblems",
)
def diffusion_map(adata, n_comps: int = 2, t: int = 1, test: bool = False):
    import numpy as np
    import scipy.sparse
    import umap

    adata = log_cp10k(adata)

    graph = umap.UMAP(transform_mode="graph").fit_transform(adata.X)
    diag_data = np.asarray(graph.sum(axis=0))
    identity = scipy.sparse.identity(graph.shape[0], dtype=np.float64)
    diag = scipy.sparse.spdiags(
        1.0 / np.sqrt(diag_data), 0, graph.shape[0], graph.shape[0]
    )
    laplacian = identity - diag * graph * diag
    num_lanczos_vectors = max(2 * n_comps + 1, int(np.sqrt(graph.shape[0])))
    eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
        laplacian,
        n_comps,
        which="SM",
        ncv=num_lanczos_vectors,
        tol=1e-4,
        v0=np.ones(laplacian.shape[0]),
        maxiter=graph.shape[0] * 5,
    )
    adata.obsm["X_emb"] = (eigenvalues**t) * eigenvectors
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
