import inspect
import os
import scprep


def r_function(filename):
    """Convert a .R file to a Python function.

    This code takes a .R file written as the body of a function.
    The function must take a single argument, ``sce``, as a
    SingleCellExperiment, and return the same object.

    Parameters
    ----------
    filename : str
        Path to the .R file relative to the calling .py file

    Returns
    -------
    fun : scprep.run.RFunction
        Python callable evaluating the R code
    """
    assert filename.endswith(".R")

    # get the path to the module that called `r_function`
    curr_frame = inspect.currentframe()
    prev_frame = inspect.getframeinfo(curr_frame.f_back)
    filepath = os.path.join(os.path.dirname(prev_frame.filename), filename)

    with open(filepath, "r") as handle:
        r_code = handle.read()

    out_fun = scprep.run.RFunction(setup="", args="sce", body=r_code)
    out_fun.__r_file__ = filepath
    return out_fun


def r_function_dual(filename):
    """Convert a .R file to a Python function.

    This code takes a .R file written as the body of a function.
    The function must take two argumens, ``sce_sc`` and ``sce_sp``, both being of type
    SingleCellExperiment, and return the same object.

    Parameters
    ----------
    filename : str
        Path to the .R file relative to the calling .py file

    Returns
    -------
    fun : scprep.run.RFunction
        Python callable evaluating the R code
    """
    assert filename.endswith(".R")

    # get the path to the module that called `r_function`
    curr_frame = inspect.currentframe()
    prev_frame = inspect.getframeinfo(curr_frame.f_back)
    filepath = os.path.join(os.path.dirname(prev_frame.filename), filename)

    with open(filepath, "r") as handle:
        r_code = handle.read()

    out_fun = scprep.run.RFunction(setup="", args=("sce_sc", "sce_sp"), body=r_code)
    out_fun.__r_file__ = filepath
    return out_fun
