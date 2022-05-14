import anndata
import numpy as np
import scipy.sparse


def split_data(
    adata: anndata.AnnData, train_frac: float = 0.9, seed: int = 0
) -> anndata.AnnData:
    """Split data using molecular cross-validation.

    Stores "train" and "test" dataset using the AnnData.obsm property.
    """
    import molecular_cross_validation.util

    random_state = np.random.RandomState(seed)

    X = adata.X

    if scipy.sparse.issparse(X):
        X = np.array(X.todense())
    if np.allclose(X, X.astype(np.int)):
        X = X.astype(np.int)
    else:
        raise TypeError("Molecular cross-validation requires integer count data.")

    X_train, X_test = molecular_cross_validation.util.split_molecules(
        X, 0.9, 0.0, random_state
    )
    adata.obsm["train"] = scipy.sparse.csr_matrix(X_train)
    adata.obsm["test"] = scipy.sparse.csr_matrix(X_test)
    return adata
