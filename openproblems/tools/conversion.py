import inspect
import os
import scprep


def r_function(filename, args="sce"):
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

    out_fun = scprep.run.RFunction(setup="", args=args, body=r_code)
    out_fun.__r_file__ = filepath
    return out_fun
