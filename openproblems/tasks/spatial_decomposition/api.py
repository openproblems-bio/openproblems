from ...data.sample import load_sample_data
from pandas.core.dtypes.common import is_categorical_dtype

import numpy as np
import pandas as pd

EPS = 1e-10


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    # test for spatial coordinates
    assert "spatial" in adata.obsm
    # check that proportions are included
    assert "proportions_true" in adata.obsm
    # make sure proportions sum to one, some precision error allowed
    assert np.all(
        np.abs(
            np.sum(
                adata[adata.obs["modality"] == "sp"].obsm["proportions_true"], axis=1
            )
            - 1  # noqa: W503
        )
        < EPS  # noqa: W503
    )
    # ensure cell type labels are found in single cell reference

    assert is_categorical_dtype(adata.obs["label"])
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "proportions_pred" in adata.obsm
    assert "proportions_true" in adata.obsm
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    # set number of cell types
    n_types = 3
    # load sample anndata
    adata = load_sample_data()
    # set spatial coordinates
    adata.obsm["spatial"] = np.random.random((adata.shape[0], 2))
    # generate proportion values
    props = pd.DataFrame(
        np.random.dirichlet(alpha=np.ones(n_types), size=adata.shape[0]),
        columns=np.arange(n_types),
        index=adata.obs.index,
    )
    adata.obsm["proportions_true"] = props
    # get anndata for single cell reference
    sc_adata = load_sample_data()
    # set labels for single cell data
    sc_adata.obs["label"] = np.random.randint(0, n_types, size=sc_adata.shape[0])
    # bind single cell reference to spatial anndata
    adata.uns["sc_reference"] = sc_adata

    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    # get number of cell types
    n_types = adata.obsm["proportions_true"].shape[1]
    # generate predicted proportions
    props = pd.DataFrame(
        np.random.dirichlet(alpha=np.ones(n_types), size=adata.shape[0]),
        columns=np.arange(n_types),
        index=adata.obs.index,
    )
    adata.obsm["proportions_pred"] = props
    return adata
