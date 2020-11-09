import numpy as np

import anndata
import warnings
import parameterized
import subprocess

from scipy import sparse


def object_name(x):
    """Get a human readable name for an object."""
    try:
        return x.__name__
    except AttributeError:
        return str(x)


def name_test(testcase_func, param_num, param):
    """Get a human readable name for a parameterized test."""
    return "%s_%s" % (
        testcase_func.__name__,
        parameterized.parameterized.to_safe_name(
            "_".join(object_name(x) for x in param.args)
        ),
    )


def ignore_warnings():
    """Ignore irrelevant warnings."""
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="is_categorical is deprecated and will be removed in a future version."
        "  Use is_categorical_dtype instead",
    )

    try:
        import numba
    except ImportError:
        return

    warnings.filterwarnings("ignore", category=numba.NumbaWarning)


def data(obsm=None):
    """Create fake data."""
    adata = anndata.AnnData(np.random.poisson(2, (100, 30)))
    if obsm is not None:
        adata.obsm[obsm] = adata.X * 2 + 1
        adata.uns["{}_obs".format(obsm)] = np.arange(adata.shape[0]) + 5
        adata.uns["{}_var".format(obsm)] = np.arange(adata.shape[1]) + 12
    return adata


def assert_array_equal(X, Y):
    """Assert two arrays to be equal, whether sparse or dense."""
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


def _format_error(process, stream):
    """Format subprocess output."""
    return "{}\nReturn code {}\n\n{}".format(
        " ".join(process.args), process.returncode, stream.decode("utf-8")
    )


def format_error_stderr(process):
    """Format subprocess output from stderr."""
    return _format_error(process, process.stderr)


def format_error_stdout(process):
    """Format subprocess output from stdout."""
    return _format_error(process, process.stdout)


def git_file_age(filename):
    """Get the age of a file's last git commit."""
    git_age = (
        run(
            ["git", "log", "-1", '--format="%ad"', "--date=unix", "--", filename],
            return_stdout=True,
        )
        .strip()
        .replace('"', "")
    )
    if git_age == "":
        return 0
    else:
        return int(git_age)


def run(
    command,
    shell=False,
    print_stdout=False,
    return_stdout=False,
    return_code=False,
    error_raises=AssertionError,
    format_error=None,
):
    """Run subprocess.

    Parameters
    ----------
    command : list of str
    shell : bool
        Run command in a new shell
    print_stdout : bool
        Print subprocess stdout to sys.stdout
    return_stdout : bool
        Return subprocess stdout
    return_code : bool
        Return subprocess exit code
    error_raises : Exception
        Which exception to raise on failure
    format_error : callable
        Function to call to generate error message. If None, chooses from
        `format_error_stderr` and `format_error_stdout` automatically.
    """
    if return_stdout and print_stdout:
        raise NotImplementedError
    elif return_stdout:
        stderr = subprocess.PIPE
        if format_error is None:
            format_error = format_error_stderr
    else:
        stderr = subprocess.STDOUT
        if format_error is None:
            format_error = format_error_stdout
    p = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=stderr)
    if print_stdout:
        while True:
            output = p.stdout.readline().decode("utf-8")
            if output == "" and p.poll() is not None:
                break
            if output:
                print(output.strip())
    p.stdout, p.stderr = p.communicate()
    output = []
    if return_stdout:
        output.append(p.stdout.decode("utf-8"))
    if return_code:
        output.append(p.returncode)
    if not return_code and not p.returncode == 0:
        raise error_raises(format_error(p))
    if output:
        return output[0] if len(output) == 1 else tuple(output)
