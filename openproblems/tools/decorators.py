from ..utils import temporary
from . import utils

import anndata
import functools
import logging
import time

log = logging.getLogger("openproblems")


def normalizer(func, *args, **kwargs):
    """Decorate a normalization function."""

    @functools.wraps(func)
    def normalize(adata, *args, obsm=None, obs=None, var=None, **kwargs):
        log.debug("Running {} normalization".format(func.__name__))
        assert isinstance(adata, anndata.AnnData)

        if obsm is not None:
            cache_name = "{}_{}".format(obsm, func.__name__)
            if cache_name in adata.obsm:
                adata.obsm[obsm] = adata.obsm[cache_name]
            else:
                obs = adata.uns[obs] if obs else adata.obs
                var = adata.uns[var] if var else adata.var
                adata_temp = anndata.AnnData(
                    adata.obsm[obsm],
                    obs=obs,
                    var=var,
                    layers={"counts": adata.obsm[obsm]},
                )
                adata_temp = func(adata_temp, *args, **kwargs)
                adata.obsm[obsm] = adata.obsm[cache_name] = adata_temp.X
        else:
            if func.__name__ in adata.layers:
                adata.X = adata.layers[func.__name__]
            else:
                adata.X = adata.layers["counts"]
                adata = func(adata, *args, **kwargs)
                adata.layers[func.__name__] = adata.X

        return adata

    return normalize


@temporary(version="1.0")
def _backport_code_version(apply_method, code_version):
    apply_method.metadata["code_version"] = code_version
    return apply_method


def method(
    method_name,
    paper_name,
    paper_reference,
    paper_year,
    code_url,
    code_version=None,
    image="openproblems",
    is_baseline=False,
):
    """Decorate a method function.

    Parameters
    ----------
    method_name : str
        Unique human readable name of the method
    paper_name : str
        Title of the seminal paper describing the method
    paper_reference : str
        BibTex key from `main.bib` referring to the paper
    paper_year : int
        Year the paper was published
    code_url : str
        Link to the code base providing the canonical implementation
    image : str, optional (default: "openproblems")
        Name of the Docker image to be used for this method
    is_baseline : bool, optional (default: False)
        If True, this method serves as a baseline for the task
    """

    def decorator(func):
        @functools.wraps(func)
        def apply_method(adata: anndata.AnnData, *args, **kwargs):
            log.debug("Running {} method".format(func.__name__))
            adata = func(adata, *args, **kwargs)
            adata.uns["is_baseline"] = is_baseline
            return adata

        apply_method.metadata = dict(
            method_name=method_name,
            paper_name=paper_name,
            paper_reference=paper_reference,
            paper_year=paper_year,
            code_url=code_url,
            image=image,
            is_baseline=is_baseline,
        )
        apply_method = _backport_code_version(apply_method, code_version)
        return apply_method

    return decorator


def metric(metric_name, maximize, paper_reference, image="openproblems"):
    """Decorate a metric function.

    Parameters
    ----------
    dataset_name : str
        Unique human readable name of the dataset

    Parameters
    ----------
    metric_name : str
        Unique human readable name of the metric
    paper_reference : str
        BibTex key from `main.bib` referring to the seminal paper in which the metric
        was defined
    maximize : bool
        If True, the metric should be maximized. If False, it should be minimized.
    image : str, optional (default: "openproblems")
        Name of the Docker image to be used for this metric
    """

    def decorator(func):
        @functools.wraps(func)
        def apply_metric(adata: anndata.AnnData, *args, **kwargs):
            log.debug("Running {} metric".format(func.__name__))
            return func(adata.copy(), *args, **kwargs)

        apply_metric.metadata = dict(
            metric_name=metric_name,
            paper_reference=paper_reference,
            maximize=maximize,
            image=image,
        )
        return apply_metric

    return decorator


def dataset(
    dataset_name=None,
    data_url=None,
    data_reference=None,
    dataset_summary=None,
    image="openproblems",
):
    """Decorate a dataset function.

    Parameters
    ----------
    dataset_name : str
        Unique human readable name of the dataset
    data_url : str
        Link to the original source of the dataset
    data_reference : str
        BibTex key from `main.bib` referring to the paper describing how the dataset was
        generated
    dataset_summary : str
        Short (<80 character) summary of the dataset
    image : str, optional (default: "openproblems")
        Name of the Docker image to be used for this dataset
    """

    def decorator(func):
        @functools.wraps(func)
        def apply_func(*args, **kwargs):
            log.debug("Loading {} dataset".format(func.__name__))
            adata = func(*args, **kwargs)
            adata.strings_to_categoricals()
            return adata

        apply_func.metadata = dict(
            dataset_name=dataset_name,
            image=image,
            data_url=data_url,
            data_reference=data_reference,
            dataset_summary=dataset_summary,
        )
        return apply_func

    return decorator


def profile(func):
    """Decorate a function for performance profiling.

    Returns
    -------
    result : dict
        Contains 'result', 'runtime_s', 'memory_mb', 'memory_leaked_mb'
    """
    import memory_profiler

    @functools.wraps(func)
    def decorated(*args, **kwargs):
        output = dict()
        utils.garbage_collection()

        def dummy():
            pass

        base_memory = memory_profiler.memory_usage(
            (dummy, tuple(), dict()),
            interval=0.1,
            max_usage=True,
        )

        def apply_func(*args, **kwargs):
            log.debug("Profiling {} function".format(func.__name__))
            start_time = time.perf_counter()
            output["result"] = func(*args, **kwargs)
            end_time = time.perf_counter()
            output["runtime_s"] = end_time - start_time

        peak_memory = memory_profiler.memory_usage(
            (apply_func, args, kwargs),
            multiprocess=True,
            include_children=True,
            interval=1,
            max_usage=True,
        )
        output["memory_mb"] = peak_memory - base_memory
        utils.garbage_collection()

        post_memory = memory_profiler.memory_usage(
            (dummy, tuple(), dict()),
            interval=0.1,
            max_usage=True,
        )
        output["memory_leaked_mb"] = max(0.0, post_memory - base_memory)
        return output

    return decorated
