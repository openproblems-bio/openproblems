import anndata
import numpy as np
import scipy.sparse


def data(sparse=False, obsm=None):
    """Create fake data."""
    data = np.random.poisson(2, (100, 30)).astype(np.float32)
    if sparse:
        data = scipy.sparse.csr_matrix(data)
    adata = anndata.AnnData(data, layers={"counts": data})
    if obsm is not None:
        adata.obsm[obsm] = adata.X * 2 + 1
        adata.uns["{}_obs".format(obsm)] = np.arange(adata.shape[0]) + 5
        adata.uns["{}_var".format(obsm)] = np.arange(adata.shape[1]) + 12
    return adata
