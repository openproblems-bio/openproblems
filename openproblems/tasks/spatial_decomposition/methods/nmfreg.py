from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp

import numpy as np


@method(
    method_name="NMF-reg",
    paper_name="Slide-seq: A scalable technology for measuring genome-wide"
    " expression at high spatial resolution",
    paper_url="https://science.sciencemag.org/content/363/6434/1463",
    paper_year=2019,
    code_url="https://github.com/tudaga/NMFreg_tutorial",
)
def nmfreg(adata, test=False, factors=None):
    """NMF-reg: NMF regression for array-based spatial transcriptomics data.

    Re-implementation from https://github.com/tudaga/NMFreg_tutorial.
    Originally developed for Slide-seq data.

    Parameters
    ----------
    adata : AnnData
        Adata with true proportions and signature matrix.

    Returns
    -------
    Adata with predicted proportions saved in `adata.obsm["proportions_pred"]`.
    """
    from scipy.optimize import nnls
    from scipy.sparse import issparse
    from sklearn.decomposition import NMF
    from sklearn.preprocessing import StandardScaler

    adata_sc, adata = split_sc_and_sp(adata)

    n_types = adata_sc.obs["label"].cat.categories.shape[0]

    factors = factors or 30

    # Learn from reference
    if issparse(adata_sc.X):
        X = adata_sc.X.toarray()
    else:
        X = adata_sc.X
    X_norm = X / X.sum(1)[:, np.newaxis]
    X_scaled = StandardScaler(with_mean=False).fit_transform(X_norm)

    model = NMF(
        n_components=factors,
        init="random",
        random_state=17,  # TODO(handle random_state)
    )
    Ha = model.fit_transform(X_scaled)
    Wa = model.components_

    cluster_df = adata.obs[["label"]].copy()
    cluster_df.loc[:, "factor"] = np.argmax(Ha, axis=1)
    cluster_df.loc[:, "code"] = cluster_df.label.values.codes
    factor_to_cluster_map = np.array(
        [
            np.histogram(
                cluster_df.loc[cluster_df.factor == k, "code"],
                bins=n_types,
                range=(0, n_types),
            )[0]
            for k in range(factors)
        ]
    ).T

    factor_to_best_celltype = np.argmax(factor_to_cluster_map, axis=0)

    factor_to_best_celltype_matrix = np.zeros((factors, n_types))
    for i, j in enumerate(factor_to_best_celltype):
        factor_to_best_celltype_matrix[i, j] = 1

    Ha_norm = StandardScaler(with_mean=False).fit_transform(Ha)
    sc_deconv = np.dot(Ha_norm**2, factor_to_best_celltype_matrix)

    sc_deconv = sc_deconv / sc_deconv.sum(1)[:, np.newaxis]

    # Start run on actual spatial data
    if issparse(adata.X):
        X_sp = adata.X.toarray()
    else:
        X_sp = adata.X
    X_sp_norm = X_sp / X_sp.sum(1)[:, np.newaxis]
    X_sp_scaled = StandardScaler(with_mean=False).fit_transform(X_sp_norm)

    bead_prop_soln = np.array(
        [nnls(Wa.T, X_sp_scaled[b, :])[0] for b in range(X_sp_scaled.shape[0])]
    )
    bead_prop_soln = StandardScaler(with_mean=False).fit_transform(bead_prop_soln)
    bead_prop = np.dot(bead_prop_soln, factor_to_best_celltype_matrix)

    prop = bead_prop / bead_prop.sum(1)[:, np.newaxis]
    adata.obsm["proportions_pred"] = prop

    adata.uns["method_code_version"] = check_version("scikit-learn")

    return adata
