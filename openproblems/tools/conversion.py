import scprep


def r_function(filename):
    """Convert a .R file to a Python function.

    This code takes a .R file written as the body of a function.
    The function must take a single argument, ``sce``, as a
    SingleCellExperiment, and return the same object.

    Parameters
    ----------
    filename : str
        Name of the .R file

    Returns
    -------
    fun : scprep.run.RFunction
        Python callable evaluating the R code
    """
    assert filename.endswith(".R")
    with open(filename, "w") as handle:
        r_code = handle.read()

    r_fun = "function evaluate(sce) {{ {code} }}".format(code=r_code)

    return scprep.run.RFunction(setup="", args="sce", body=r_fun)
