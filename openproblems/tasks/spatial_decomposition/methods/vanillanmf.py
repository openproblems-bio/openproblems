from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_version

from sklearn.decomposition import NMF
import numpy as np
import pandas as pd


@method(
    method_name="Non-Negative Matrix Factorization (NMF).",
    paper_name="Fast local algorithms for large scale nonnegative matrix and tensor factorizations", 
    paper_url="https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.214.6398&rep=rep1&type=pdf",
    paper_year=2009,
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html",
    code_version=check_version("scikit-learn"),
)

# Used from https://github.com/theislab/scanpy/issues/181#issuecomment-534867254
def grouped_obs_mean(adata, group_key, layer=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out.loc[:,group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out.T


def nmf_raw(adata):
    n_types = adata.obsm["proportions_true"].size[1]
    adata_sc = adata.uns["sc_reference"].copy()

    vanila_nmf_model = NMF(
        n_components=n_types,
        beta_loss='kullback-leibler', 
        solver='mu', 
        max_iter=4000, 
        alpha=.1,
        init='custom', 
        random_state=17
    )

    # Make profiles from single-cell expression 
    # dataset
    profile_mean = grouped_obs_mean(adata_sc, 'celltype')

    Wa = vanila_nmf_model.fit_transform(
        H=profile_mean.values
    )

    prop = Wa_norm = Wa / Wa.sum(1)[:,np.newaxis]
    adata.obsm['proportions_pred'] = prop
    
    return adata