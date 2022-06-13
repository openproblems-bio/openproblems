import warnings


def ignore_warnings():
    """Ignore irrelevant warnings."""
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="is_categorical is deprecated and will be removed in a future version."
        "  Use is_categorical_dtype instead",
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
        message="The global conversion available with activate() is deprecated",
        module="rpy2.robjects.numpy2ri",
    )
    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message="The global conversion available with activate() is deprecated",
        module="rpy2.robjects.pandas2ri",
    )

    try:
        import numba
    except ImportError:
        return

    warnings.filterwarnings("ignore", category=numba.NumbaWarning)
