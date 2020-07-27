from decorator import decorator
import anndata


@decorator
def normalizer(func, *args, **kwargs):
    adata = args[0]
    assert isinstance(adata, anndata.AnnData)
    if "store" in kwargs:
        store = kwargs["store"]
        del kwargs["store"]
    else:
        store = False

    if func.__name__ in adata.layers:
        adata.X = adata.layers[func.__name__]
    else:
        func(*args, **kwargs)
        if store:
            adata.layers[func.__name__] = adata.X
