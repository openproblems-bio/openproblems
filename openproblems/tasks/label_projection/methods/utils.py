def pca_op(adata_train, adata_test, n_components=100):
    import scipy.sparse
    import sklearn.decomposition

    is_sparse = scipy.sparse.issparse(adata_train.X)

    min_components = min(
        [adata_train.shape[0], adata_test.shape[0], adata_train.shape[1]]
    )
    if is_sparse:
        min_components -= 1
    n_components = min([n_components, min_components])
    if is_sparse:
        pca_op = sklearn.decomposition.TruncatedSVD
    else:
        pca_op = sklearn.decomposition.PCA
    return pca_op(n_components=n_components)
