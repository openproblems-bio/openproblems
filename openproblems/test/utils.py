import numpy as np

import openproblems
import anndata
import warnings
import parameterized
import subprocess
import time
import threading
import queue
import sys

import scipy.sparse
import packaging.version


def object_name(x):
    """Get a human readable name for an object."""
    if hasattr(x, "__name__"):
        return x.__name__
    elif hasattr(x, "__func__"):
        return object_name(x.__func__)
    else:
        return str(x)


def name_test(testcase_func, param_num, param):
    """Get a human readable name for a parameterized test."""
    args = param.values() if isinstance(param, dict) else param.args
    return "%s_%s" % (
        testcase_func.__name__,
        parameterized.parameterized.to_safe_name(
            "_".join(object_name(x) for x in args)
        ),
    )


def assert_finite(X):
    """Assert numpy or scipy matrix is finite."""
    X = X.data if scipy.sparse.issparse(X) else X
    assert np.all(np.isfinite(X))
    return True


def future_warning(msg, error_version, error_category, warning_category=FutureWarning):
    """Raise a warning until a specific version, then raise an error."""
    current_version = packaging.version.parse(openproblems.__version__)
    if current_version < packaging.version.parse(error_version):
        msg += " This will raise a {} in openproblems v{}".format(
            error_category.__name__, error_version
        )
        warnings.warn(msg, warning_category)
    else:
        raise error_category(msg)


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
    if scipy.sparse.issparse(X) and scipy.sparse.issparse(Y):
        X = X.tocsr()
        Y = Y.tocsr()
        np.testing.assert_array_equal(X.data, Y.data)
        np.testing.assert_array_equal(X.indices, Y.indices)
        np.testing.assert_array_equal(X.indptr, Y.indptr)
    else:
        X = np.asarray(X)
        Y = np.asarray(Y)
        np.testing.assert_array_equal(X, Y)


def format_error_timeout(process, timeout, stream):
    """Format subprocess output on timeout."""
    return "{}\nTimed out after {} s\n\n{}".format(
        " ".join(process.args),
        timeout,
        NonBlockingStreamReader(stream).read().decode("utf-8"),
    )


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
    timeout=3600,
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
    if timeout is not None:
        runtime = 0
        if p.poll() is None:
            time.sleep(1)
            runtime += 1
        if runtime > timeout:
            raise RuntimeError(
                format_error_timeout(
                    p, timeout, p.stderr if stderr is subprocess.PIPE else p.stdout
                )
            )

    if print_stdout:
        while True:
            output = p.stdout.readline().decode("utf-8")
            if output == "" and p.poll() is not None:
                break
            if output:
                print(output.strip())
                sys.stdout.flush()
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


class NonBlockingStreamReader:
    """Filestream that can be read without blocking.

    Parameters
    ----------
    stream: the stream to read from.
        Usually a process' stdout or stderr.
    """

    @staticmethod
    def populateQueue(stream, queue):
        """Collect lines from 'stream' and put them in 'queue'."""
        while True:
            line = stream.readline()
            if line:
                queue.put(line)
            else:
                raise RuntimeError("Stream ended without EOF.")

    def __init__(self, stream):
        """Initialize class."""
        self.stream = stream
        self.queue = queue.Queue()

        self.thread = threading.Thread(
            target=self.populateQueue, args=(self.stream, self.queue)
        )
        self.thread.daemon = True
        self.thread.start()  # start collecting lines from the stream

    def readline(self, timeout=None):
        """Read one line from stream if available.

        Parameters
        ----------
        timeout : int, optional (default: None)
            Seconds to wait for output. If None, block until response.
        """
        try:
            return self.queue.get(block=timeout is not None, timeout=timeout)
        except queue.Empty:
            return None

    def read(self, timeout=None):
        """Read all available lines from stream.

        Parameters
        ----------
        timeout : int, optional (default: None)
            Seconds to wait for output. If None, block until EOF.
        """
        output = b""
        while True:
            next_line = self.readline(timeout=timeout)
            if next_line is None:
                return output
            output += next_line + b"\n"
