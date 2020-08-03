import warnings
import parameterized
import numpy as np
import anndata


def object_name(x):
    try:
        return x.__name__
    except AttributeError:
        return str(x)


def name_test(testcase_func, param_num, param):
    return "%s_%s" % (
        testcase_func.__name__,
        parameterized.parameterized.to_safe_name(
            "_".join(object_name(x) for x in param.args)
        ),
    )


def ignore_numba_warnings():
    try:
        import numba
    except ImportError:
        return

    warnings.filterwarnings("ignore", category=numba.NumbaWarning)


def data(obsm=None):
    adata = anndata.AnnData(np.random.poisson(2, (100, 10)))
    if obsm is not None:
        adata.obsm[obsm] = adata.X * 2 + 1
        adata.uns["{}_obs".format(obsm)] = np.arange(adata.shape[0]) + 5
        adata.uns["{}_var".format(obsm)] = np.arange(adata.shape[1]) + 12
    return adata
