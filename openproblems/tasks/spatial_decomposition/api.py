from ...data.sample import load_sample_data
from ._utils import merge_sc_and_sp
from pandas.core.dtypes.common import is_categorical_dtype

import numpy as np

EPS = 1e-10


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    # test for spatial coordinates
    assert "spatial" in adata.obsm
    assert isinstance(adata.obsm["spatial"], np.ndarray)
    # check that proportions are included
    assert "proportions_true" in adata.obsm
    assert isinstance(adata.obsm["proportions_true"], np.ndarray)
    # make sure proportions sum to one, some precision error allowed
    proportions_sum = np.sum(
        adata[adata.obs["modality"] == "sp"].obsm["proportions_true"], axis=1
    )
    np.testing.assert_allclose(proportions_sum, 1)
    # ensure cell type labels are found in single cell reference
    assert is_categorical_dtype(adata.obs["label"])
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "proportions_pred" in adata.obsm
    assert isinstance(adata.obsm["proportions_pred"], np.ndarray)
    assert "proportions_true" in adata.obsm
    assert isinstance(adata.obsm["proportions_true"], np.ndarray)
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    # set number of cell types
    n_types = 3
    # load sample anndata
    adata_spatial = load_sample_data()
    # modify index
    adata_spatial.obs.index = "spatial_" + adata_spatial.obs.index
    # set spatial coordinates
    adata_spatial.obsm["spatial"] = np.random.random((adata_spatial.shape[0], 2))
    # generate proportion values
    props = np.random.dirichlet(alpha=np.ones(n_types), size=adata_spatial.shape[0])
    adata_spatial.obsm["proportions_true"] = props
    # get anndata for single cell reference
    adata_sc = load_sample_data()
    # set labels for single cell data
    adata_sc.obs["label"] = np.random.randint(
        0, n_types, size=adata_sc.shape[0]
    ).astype(str)
    # bind single cell reference to spatial anndata
    adata_merged = merge_sc_and_sp(adata_sc, adata_spatial)
    return adata_merged


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    # get number of cell types
    n_types = adata.obsm["proportions_true"].shape[1]
    # generate predicted proportions
    props = np.random.dirichlet(alpha=np.ones(n_types), size=adata.shape[0])
    adata.obsm["proportions_pred"] = props
    return adata
