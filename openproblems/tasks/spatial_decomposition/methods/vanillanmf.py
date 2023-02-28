from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import obs_means
from ..utils import split_sc_and_sp


@method(
    method_name="Non-Negative Matrix Factorization (NMF)",
    paper_name=(
        "Fast local algorithms for large scale nonnegative "
        "matrix and tensor factorizations"
    ),
    paper_reference="cichocki2009fast",
    paper_year=2009,
    code_url=(
        "https://scikit-learn.org/stable/modules/generated/"
        "sklearn.decomposition.NMF.html"
    ),
)
def nmf(adata, test=False, max_iter=None, random_state=17):
    """NMF for spatial deconvolution."""
    from scipy.sparse import issparse
    from sklearn.decomposition import NMF

    import numpy as np

    adata_sc, adata = split_sc_and_sp(adata)
    n_types = adata_sc.obs["label"].cat.categories.shape[0]

    if test:
        max_iter = max_iter or 10
    else:  # pragma: nocover
        max_iter = max_iter or 4000

    vanila_nmf_model = NMF(
        n_components=n_types,
        beta_loss="kullback-leibler",
        solver="mu",
        max_iter=max_iter,
        alpha_W=0.1,
        alpha_H=0.1,
        init="custom",
        random_state=random_state,
    )

    # Make profiles from single-cell expression dataset
    adata_means = obs_means(adata_sc, "label")

    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    Wa = vanila_nmf_model.fit_transform(
        X.astype(adata_means.X.dtype),
        H=adata_means.X,
        W=np.ones((adata.shape[0], n_types), dtype=adata_means.X.dtype),
    )

    prop = Wa / Wa.sum(1)[:, np.newaxis]
    adata.obsm["proportions_pred"] = prop

    adata.uns["method_code_version"] = check_version("scikit-learn")

    return adata
