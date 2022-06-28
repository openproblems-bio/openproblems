from ....tools.decorators import method
from ....tools.utils import check_version

# import numpy as np
import scanpy as sc


def remove_zeros(Xmat):
    if scipy.sparse.issparse(Xmat):
        Xmat = Xmat.tocsc()
        # creates boolean list of missing
        missing_bool = np.apply_along_axis(np.count_nonzero, 0, Xmat) == 0
        # removes missing genes
        Xmat[, (np.apply_along_axis(np.count_nonzero, 0, Xmat) == 0)]
        # creates list of true (missing) indices
        inds = list(compress(xrange(len(missing_bool)), missing_bool))
    return Xmat, inds

# from https://stackoverflow.com/questions/53870310


def insert_at(arr, output_size, indices):
    result = np.zeros(output_size)
    existing_indices = [np.setdiff1d(np.arange(axis_size), 
                                     axis_indices, assume_unique=True)
                        for axis_size, axis_indices in zip(output_size, indices)]
    result[np.ix_(*existing_indices)] = arr
    return result


def _dca(adata, test=False, epochs=None):
    if test:
        epochs = 30
    else:
        epochs = epochs or 300
    from dca.api import dca

    # remove zeros from train
    Xmat, inds = remove_zeros(adata.obsm["train"])
    # make adata object with train counts
    adata2 = sc.AnnData(Xmat)
    # run DCA
    dca(adata2, epochs=epochs)
    Xmat = adata2.X
    # insert zero rows back into Xmat at inds
    Xmat = insert_at(Xmat,
                     (adata.obsm["train"].shape[0], adata.obsm["train"].shape[1]),
                     (0, list(inds)))[1:, :]
    # set denoised to Xmat
    adata.obsm["denoised"] = Xmat
    # check version of dca
    adata.uns["method_code_version"] = check_version("dca")
    return adata


@method(
    method_name="DCA",
    paper_name="Single-cell RNA-seq denoising using...",
    paper_url="https://www.nature.com/articles/s41467-018-07931-2",
    paper_year=2019,
    code_url="https://github.com/theislab/dca",
    image="openproblems-python-tf2.4",
)
def dca(adata, test=False, epochs=None):
    return _dca(adata, test=test, epochs=epochs)
