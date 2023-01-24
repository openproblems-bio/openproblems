import anndata
import functools
import scanpy as sc
import scprep

def normalizer(func, *args, **kwargs):
    """Decorate a normalization function."""

    @functools.wraps(func)
    def normalize(adata, *args, obsm=None, obs=None, var=None, **kwargs):
        # log.debug("Running {} normalization".format(func.__name__))
        assert isinstance(adata, anndata.AnnData)

        if obsm is not None:
            cache_name = "{}_{}".format(obsm, func.__name__)
            if cache_name in adata.obsm:
                adata.obsm[obsm] = adata.obsm[cache_name]
            else:
                obs = adata.uns[obs] if obs else adata.obs
                var = adata.uns[var] if var else adata.var
                adata_temp = anndata.AnnData(adata.obsm[obsm], obs=obs, var=var)
                func(adata_temp, *args, **kwargs)
                adata.obsm[obsm] = adata.obsm[cache_name] = adata_temp.X
        else:
            if func.__name__ in adata.layers:
                adata.X = adata.layers[func.__name__]
            else:
                func(adata, *args, **kwargs)
                adata.layers[func.__name__] = adata.X

    return normalize

def _cpm(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")

@normalizer
def cpm(adata):
    """Normalize data to counts per million."""
    _cpm(adata)

@normalizer
def log_cpm(adata):
    """Normalize data to log counts per million."""
    _cpm(adata)
    sc.pp.log1p(adata)

@normalizer
def sqrt_cpm(adata):
    """Normalize data to sqrt counts per million."""
    _cpm(adata)
    adata.X = scprep.transform.sqrt(adata.X)
