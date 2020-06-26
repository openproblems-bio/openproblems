import mnnpy
import numpy as np
from sklearn.decomposition import TruncatedSVD


def mnn(adata, n_svd=100):
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = TruncatedSVD(n_svd).fit_transform(adata.uns["mode2"].X)
    (X_corrected, Y_corrected), _, _ = mnnpy.mnn_correct(
        X_pca, Y_pca, var_index=np.arange(n_svd), do_concatenate=False
    )
    adata.obsm["aligned"] = X_corrected
    adata.uns["mode2"].obsm["aligned"] = Y_corrected
