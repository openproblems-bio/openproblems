from decorator import decorator
import anndata


@decorator
def normalizer(func, *args, **kwargs):
    adata = args[0]
    assert isinstance(adata, anndata.AnnData)

    if func.__name__ in adata.layers:
        adata.X = adata.layers[func.__name__]
    else:
        func(*args, **kwargs)
        adata.layers[func.__name__] = adata.X
