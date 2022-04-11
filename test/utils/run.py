from . import streams

import logging
import subprocess
import sys
import time

log = logging.getLogger("openproblems")


def format_error_timeout(process, timeout, stream):
    """Format subprocess output on timeout."""
    return "{}\nTimed out after {} s\n\n{}".format(
        " ".join(process.args),
        timeout,
        streams.NonBlockingStreamReader(stream).read().decode("utf-8"),
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


def _run_failed(process, error_raises, format_error):
    raise error_raises(format_error(process))


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

    log.debug("Running subprocess: {}".format(command))
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

    log.debug("Awaiting subprocess completion")
    if print_stdout:
        while True:
            output = p.stdout.readline().decode("utf-8")
            if output == "" and p.poll() is not None:
                break
            if output:
                print(output.strip())
                sys.stdout.flush()
    else:
        p.wait()

    log.debug("Subprocess complete")
    p.stdout, p.stderr = p.communicate()
    output = []
    if return_stdout:
        output.append(p.stdout.decode("utf-8"))
    if return_code:
        output.append(p.returncode)
    if not return_code and not p.returncode == 0:
        _run_failed(p, error_raises, format_error)
    if output:
        return output[0] if len(output) == 1 else tuple(output)
