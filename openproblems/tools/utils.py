from decorator import decorator
import anndata


@decorator
def normalizer(func, *args, obsm=None, obs=None, var=None, **kwargs):
    adata = args[0]
    assert isinstance(adata, anndata.AnnData)

    if obsm is not None:
        cache_name = "{}_{}".format(obsm, func.__name__)
        if cache_name in adata.obsm:
            adata.obsm[obsm] = adata.obsm[cache_name]
        else:
            obs = obs if obs else adata.obs
            var = var if var else adata.var
            adata_temp = anndata.AnnData(adata.obsm[obsm], obs=obs, var=var)
            func(adata_temp, *args[1:], **kwargs)
            adata.obsm[obsm] = adata.obsm[cache_name] = adata_temp.X
    else:
        if func.__name__ in adata.layers:
            adata.X = adata.layers[func.__name__]
        else:
            func(*args, **kwargs)
            adata.layers[func.__name__] = adata.X
