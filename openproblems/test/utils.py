import numpy as np

import anndata
import warnings
import parameterized
import subprocess

from scipy import sparse


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


def ignore_warnings():
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="is_categorical is deprecated and will be removed in a future version.  Use is_categorical_dtype instead",
    )

    try:
        import numba
    except ImportError:
        return

    warnings.filterwarnings("ignore", category=numba.NumbaWarning)


def data(obsm=None):
    adata = anndata.AnnData(np.random.poisson(2, (100, 30)))
    if obsm is not None:
        adata.obsm[obsm] = adata.X * 2 + 1
        adata.uns["{}_obs".format(obsm)] = np.arange(adata.shape[0]) + 5
        adata.uns["{}_var".format(obsm)] = np.arange(adata.shape[1]) + 12
    return adata


def assert_array_equal(X, Y):
    assert X.shape == Y.shape
    if sparse.issparse(X) and sparse.issparse(Y):
        X = X.tocsr()
        Y = Y.tocsr()
        np.testing.assert_array_equal(X.data, Y.data)
        np.testing.assert_array_equal(X.indices, Y.indices)
        np.testing.assert_array_equal(X.indptr, Y.indptr)
    else:
        X = np.asarray(X)
        Y = np.asarray(Y)
        np.testing.assert_array_equal(X, Y)


def format_error_stderr(process):
    return "Return code {}\n\n{}".format(
        process.returncode, process.stdout.decode("utf-8")
    )


def format_error_stdout(process):
    return "Return code {}\n\n{}".format(
        process.returncode, process.stdout.decode("utf-8")
    )


def run(
    command,
    shell=False,
    return_stdout=False,
    return_code=False,
    error_raises=AssertionError,
    format_error=None,
):
    if return_stdout:
        stderr = subprocess.PIPE
        if format_error is None:
            format_error = format_error_stderr
    else:
        stderr = subprocess.STDOUT
        if format_error is None:
            format_error = format_error_stdout
    p = subprocess.run(command, shell=shell, stdout=subprocess.PIPE, stderr=stderr)
    output = []
    if return_stdout:
        output.append(p.stdout.decode("utf-8"))
    if return_code:
        output.append(p.returncode)
    if not return_code and not p.returncode == 0:
        raise error_raises(format_error(p))
    if output:
        return tuple(output)
