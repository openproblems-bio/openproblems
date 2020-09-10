import functools
from memory_profiler import memory_usage
import time
import gc


def method(method_name, paper_name, paper_url, paper_year, code_url):
    def decorator(func):
        @functools.wraps(func)
        def apply_method(*args, **kwargs):
            return func(*args, **kwargs)

        apply_method.metadata = dict(
            method_name=method_name,
            paper_name=paper_name,
            paper_url=paper_url,
            paper_year=paper_year,
            code_url=code_url,
        )
        return apply_method

    return decorator


def metric(metric_name, maximize):
    def decorator(func):
        @functools.wraps(func)
        def apply_metric(*args, **kwargs):
            return func(*args, **kwargs)

        apply_metric.metadata = dict(metric_name=metric_name, maximize=maximize)
        return apply_metric

    return decorator


def profile(func):
    @functools.wraps(func)
    def decorated(*args, **kwargs):
        output = dict()
        gc.collect()

        def dummy():
            pass

        base_memory = memory_usage(
            (dummy, tuple(), dict()),
            interval=0.1,
            max_usage=True,
        )

        def apply_func(*args, **kwargs):
            start_time = time.perf_counter()
            output["result"] = func(*args, **kwargs)
            end_time = time.perf_counter()
            output["runtime_s"] = end_time - start_time

        peak_memory = memory_usage(
            (apply_func, args, kwargs),
            multiprocess=True,
            include_children=True,
            interval=1,
            max_usage=True,
        )
        output["memory_mb"] = peak_memory - base_memory
        gc.collect()

        post_memory = memory_usage(
            (dummy, tuple(), dict()),
            interval=0.1,
            max_usage=True,
        )
        output["memory_leaked_mb"] = post_memory - base_memory
        return output

    return decorated
