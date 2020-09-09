import functools


def method(method_name, paper_name, paper_url, paper_year, code_url):
    def decorator(func, *args, **kwargs):
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
    def decorator(func, *args, **kwargs):
        @functools.wraps(func)
        def apply_metric(*args, **kwargs):
            return func(*args, **kwargs)

        apply_metric.metadata = dict(metric_name=metric_name, maximize=maximize)
        return apply_metric

    return decorator
