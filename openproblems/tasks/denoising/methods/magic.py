from ....tools.decorators import baseline_method
from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep

_magic_method = functools.partial(
    method,
    method_summary=(
        "MAGIC (Markov Affinity-based Graph Imputation of Cells) is a method for"
        " imputation and denoising of noisy or dropout-prone single cell RNA-sequencing"
        " data. Given a normalised scRNA-seq expression matrix, it first calculates"
        " Euclidean distances between each pair of cells in the dataset, which is then"
        " augmented using a Gaussian kernel (function) and row-normalised to give a"
        " normalised affinity matrix. A t-step markov process is then calculated, by"
        " powering this affinity matrix t times. Finally, the powered affinity matrix"
        " is right-multiplied by the normalised data, causing the final imputed values"
        " to take the value of a per-gene average weighted by the affinities of cells."
        " The resultant imputed matrix is then rescaled, to more closely match the"
        " magnitude of measurements in the normalised (input) matrix."
    ),
    paper_name=(
        "Recovering Gene Interactions from Single-Cell Data Using Data Diffusion"
    ),
    paper_reference="van2018recovering",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    image="openproblems-python-extras",
)


def _magic(adata, solver, normtype="sqrt", reverse_norm_order=False, **kwargs):
    from magic import MAGIC

    if normtype == "sqrt":
        norm_fn = np.sqrt
        denorm_fn = np.square
    elif normtype == "log":
        norm_fn = np.log1p
        denorm_fn = np.expm1
    else:
        raise NotImplementedError

    X = adata.obsm["train"]
    if reverse_norm_order:
        # inexplicably, this sometimes performs better
        X = scprep.utils.matrix_transform(X, norm_fn)
        X, libsize = scprep.normalize.library_size_normalize(
            X, rescale=1, return_library_size=True
        )
    else:
        X, libsize = scprep.normalize.library_size_normalize(
            X, rescale=1, return_library_size=True
        )
        X = scprep.utils.matrix_transform(X, norm_fn)

    Y = MAGIC(solver=solver, **kwargs, verbose=False).fit_transform(
        X, genes="all_genes"
    )

    Y = scprep.utils.matrix_transform(Y, denorm_fn)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y

    adata.uns["method_code_version"] = check_version("magic-impute")
    return adata


@_magic_method(
    method_name="MAGIC",
)
def magic(adata, test=False):
    return _magic(adata, solver="exact", normtype="sqrt")


@_magic_method(
    method_name="MAGIC (reversed normalization)",
)
def magic_reverse_norm(adata, test=False):
    return _magic(adata, solver="exact", normtype="sqrt", reverse_norm_order=True)


@_magic_method(
    method_name="MAGIC (approximate)",
)
def magic_approx(adata, test=False):
    return _magic(adata, solver="approximate", normtype="sqrt")


@_magic_method(
    method_name="MAGIC (approximate, reversed normalization)",
)
def magic_approx_reverse_norm(adata, test=False):
    return _magic(adata, solver="approximate", normtype="sqrt", reverse_norm_order=True)


@baseline_method(
    method_name="KNN smoothing",
    method_summary=(
        "KNN-smoothing is a method for denoising data based on the k-nearest"
        " neighbours. Given a normalised scRNA-seq matrix, KNN-smoothing calculates a"
        " k-nearest neighbour matrix using Euclidean distances between cell pairs. Each"
        " cell’s denoised expression is then defined as the average expression of each"
        " of its neighbours."
    ),
    is_baseline=False,
    image="openproblems-python-extras",
)
def knn_naive(adata, test=False):
    return _magic(adata, solver="exact", normtype="log", decay=None, t=1)
