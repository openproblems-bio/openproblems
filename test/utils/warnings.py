import warnings


def ignore_warnings():
    """Ignore irrelevant warnings."""
    warnings.simplefilter("error")
    warnings.filterwarnings(
        "always",
        category=FutureWarning,
        message="self.adata.X is loaded as dense",
    )
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="The `inplace` parameter in pandas.Categorical.reorder_categories"
        " is deprecated",
        module="anndata._core.anndata",
    )

    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message="This package has been superseded by the `leidenalg` package",
        module="scanpy.tools._louvain",
    )
    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message="Please use `spmatrix` from the `scipy.sparse` namespace",
        module="anndata._core.merge",
    )
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="The `inplace` parameter in pandas.Categorical.reorder_categories",
        module="anndata._core.anndata",
    )
    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message=r"The global conversion available with activate\(\) is deprecated",
        module="rpy2.robjects.numpy2ri",
    )
    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message=r"The global conversion available with activate\(\) is deprecated",
        module="rpy2.robjects.pandas2ri",
    )

    try:
        import numba

        warnings.filterwarnings("ignore", category=numba.NumbaWarning)
    except ImportError:
        pass
