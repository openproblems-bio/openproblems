import anndata
import pandas as pd
import numpy as np


def generate_fake_anndata():
    X = np.zeros((3, 3))
    obs = pd.DataFrame([['celseq', 'gamma', 0.02], ['celseq', 'gamma', 0.07], ['celseq', 'gamma', 0.03]],
                       columns=['tech', 'celltype', 'size_factors'],
                       index=['A1BG', 'A1CF', 'A2M'])
    var = pd.DataFrame(index=['A1BG', 'A1CF', 'A2M'])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.layers['counts'] = X
    return adata
