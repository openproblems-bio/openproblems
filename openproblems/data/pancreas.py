import os
import tempfile

import numpy as np
import scprep
import scanpy as sc

from .utils import loader
from ..tools import normalize as n

URL = (
    "https://ndownloader.figshare.com/files/22891151"
)


@loader
def load_pancreas(test=False, preprocess=True):
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False)

        # Subsample pancreas data
        cts = ['delta', 'gamma']
        batches = ['inDrop2', 'smarter', 'celseq']
        adata = adata[(adata.obs['batch'].isin(batches)) &
                      (adata.obs['labels'].isin(cts)), :200].copy()
        # Note: could also use 200-500 HVGs rather than 200 random genes

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "pancreas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Correct covariate naming and remove preprocessing
            prep_pancreas(adata)

            #Preprocessing
            if preprocess:
                for fn in dir(n):
                    adata.layers[str(fn)] = fn(adata)
                
        return adata



def prep_pancreas(adata):
    # Rename categories
    adata.obs['batch'] = adata.obs['tech'].copy()
    adata.obs['labels'] = adata.obs['celltype'].copy()
    
    # Assign training/test
    test_batches = ['inDrop4', 'smartseq2']
    adata.obs['is_train'] = [False if adata.obs['batch'][idx] in test_batches
                             else True
                             for idx in adata.obs_names]

    # Drop excess covariates
    adata.obs = adata.obs.drop(columns=['tech', 'size_factors', 'celltype'])

    # Remove processing
    adata.X = adata.layers['counts']
    del adata.layers['counts']
