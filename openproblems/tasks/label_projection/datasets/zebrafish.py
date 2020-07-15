import os

import numpy as np
import anndata
import scprep
import scanpy as sc

def zebrafish(test=False, split='lab'):
    if split not in ['lab', 'random']:
        raise NotImplementedError("split must be in ['lab', 'random']")

    url = 'https://ndownloader.figshare.com/files/23820602?private_link=a23de5c2117a69d64af4'
    home = os.path.expanduser('~/')
    destination = os.path.join(home, 'zebrafish.h5ad')
    if not os.path.exists(destination):
        scprep.io.download.download_url(url, destination)
    adata = anndata.read_h5ad(destination)
    adata.obs["labels"] = adata.obs['cell_type']

    if split == 'lab':
        adata.obs["is_train"] = adata.obs['lab'] == 'Schier'
    elif split == 'random':
        adata.obs["is_train"] = np.random.choice(
            [True, False], adata.shape[0], replace=True
        )

    if test:
        sc.pp.subsample(adata, n_obs=100)
    return adata

def zebrafish_lab(test=False):
    return zebrafish(test, 'lab')

def zebrafish_random(test=False):
    return zebrafish(test, 'random')
