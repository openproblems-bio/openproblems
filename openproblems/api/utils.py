import openproblems
import os


class NoSuchFunctionError(RuntimeError):
    pass


def module_to_str(module):
    """Convert a Python module to a task name."""
    return module.__name__.split(".")[-1]


def str_to_task(task_name):
    """Convert a task name to a module."""
    return getattr(openproblems.tasks, task_name)


def function_to_str(fun):
    """Convert a function to a function name."""
    return fun.__name__


def get_function(task_name, function_type, function_name):
    """Get a function from a task."""
    # check function type
    function_type = function_type.lower()
    assert function_type in ["datasets", "methods", "metrics", "api"]

    # get function
    task = str_to_task(task_name)
    functions = getattr(task, function_type)
    try:
        fun = getattr(functions, function_name)
        assert callable(fun)
    except (AssertionError, AttributeError):
        raise NoSuchFunctionError(
            "Task {} has no {} '{}'".format(
                task_name, function_type[:-1], function_name
            )
        )
    return fun


def print_output(output):
    """Print output of the CLI."""
    if output is None:
        pass
    elif isinstance(output, list):
        print("\n".join(output))
    else:
        print(output)


def write_h5ad(adata, filename):
    if os.path.isfile(filename):
        os.unlink(filename)
    adata.write_h5ad(filename)
