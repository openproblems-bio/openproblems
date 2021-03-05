import os
import scprep


def r_function(filename, caller):
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
    assert filename.endswith(".R")
    filepath = os.path.join(os.path.dirname(caller), filename)
    with open(filepath, "w") as handle:
        r_code = handle.read()

    r_fun = "function evaluate(sce) {{ {code} }}".format(code=r_code)

    return scprep.run.RFunction(setup="", args="sce", body=r_fun)
