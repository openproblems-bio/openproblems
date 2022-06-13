from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import obs_means
from .._utils import split_sc_and_sp


@method(
    method_name="Non-Negative Matrix Factorization (NMF).",
    paper_name="Fast local algorithms for large scale nonnegative matrix and tensor factorizations",  # noqa: E501
    paper_url="https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.214.6398&rep=rep1&type=pdf",  # noqa: E501
    paper_year=2009,
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html",  # noqa: E501
    code_version=check_version("scikit-learn"),
)
def nmf(adata, test=False):
    """NMF for spatial deconvolution."""
    from scipy.sparse import issparse
    from sklearn.decomposition import NMF

    import numpy as np

    adata_sc, adata = split_sc_and_sp(adata)
    n_types = adata_sc.obs["label"].cat.categories.shape[0]

    vanila_nmf_model = NMF(
        n_components=n_types,
        beta_loss="kullback-leibler",
        solver="mu",
        max_iter=4000,
        alpha=0.1,
        init="custom",
        random_state=17,  # TODO(handle random_state)
    )

    # Make profiles from single-cell expression dataset
    adata_means = obs_means(adata_sc, "label")

    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    Wa = vanila_nmf_model.fit_transform(
        X, H=adata_means.X, W=np.ones((adata.shape[0], n_types), dtype=np.float32)
    )

    prop = Wa / Wa.sum(1)[:, np.newaxis]
    adata.obsm["proportions_pred"] = prop

    return adata
