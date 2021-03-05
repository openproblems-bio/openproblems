import os
import scprep
import inspect


def r_function(filename):
    """Convert a .R file to a Python function.

    This code takes a .R file written as the body of a function.
    The function must take a single argument, ``sce``, as a
    SingleCellExperiment, and return the same object.

    Parameters
    ----------
    filename : str
        Name of the .R file
    caller : str
        Path to the calling .py file (obtained with ``__file__``)

    Returns
    -------
    fun : scprep.run.RFunction
        Python callable evaluating the R code
    """
    curr_frame = inspect.currentframe()
    prev_frame = inspect.getframeinfo(curr_frame.f_back)
    assert filename.endswith(".R")
    filepath = os.path.join(os.path.dirname(prev_frame.filename), filename)
    with open(filepath, "r") as handle:
        r_code = handle.read()

    r_fun = "function evaluate(sce) {{ {code} }}".format(code=r_code)

    return scprep.run.RFunction(setup="", args="sce", body=r_fun)
