from . import TEMPDIR

import parameterized


def object_name(x):
    """Get a human readable name for an object."""
    if hasattr(x, "__name__"):
        return x.__name__
    elif hasattr(x, "__func__"):
        return object_name(x.__func__)
    else:
        return str(x)


def name_test(testcase_func, param_num, param):
    """Get a human readable name for a parameterized test."""
    args = param.values() if isinstance(param, dict) else param.args

    return "%s_%s" % (
        testcase_func.__name__,
        parameterized.parameterized.to_safe_name(
            "_".join(object_name(x) for x in args if x != TEMPDIR.name)
        ),
    )
