from . import utils


def list_functions(task_name, function_type):
    """Get the name of the functions associated with a task."""
    # check function type
    function_type = function_type.upper()
    assert function_type in ["DATASETS", "METHODS", "METRICS"]

    # print functions
    task = utils.str_to_task(task_name)
    functions = [utils.function_to_str(fun) for fun in getattr(task, function_type)]
    return functions


def main(args):
    """Run the ``list`` subcommand."""
    functions = list_functions(args.task, args.function_type)
    return functions
