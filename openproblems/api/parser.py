import argparse
import pathlib


def filepath(string):
    """Require an argument type to be a valid path."""
    if string is None:
        return string
    path = pathlib.Path(string)
    if path.exists and path.is_dir():
        raise argparse.ArgumentTypeError("Value must not be a directory")
    if not path.parent.is_dir():
        try:
            path.parent.mkdir(exist_ok=True, parents=True)
        except OSError as e:
            raise argparse.ArgumentTypeError(
                "Failed to create directory {}: {}".format(path.parent, e)
            )
    return path


def parse_function_type(parser):
    """Create a mutually exclusive selection group for datasets, methods and metrics."""
    selection_group = parser.add_mutually_exclusive_group(required=True)
    selection_group.add_argument(
        "--datasets",
        dest="function_type",
        action="store_const",
        const="DATASETS",
        help="List datasets from the selected task",
    )
    selection_group.add_argument(
        "--methods",
        dest="function_type",
        action="store_const",
        const="METHODS",
        help="List methods from the selected task",
    )
    selection_group.add_argument(
        "--metrics",
        dest="function_type",
        action="store_const",
        const="METRICS",
        help="List metrics from the selected task",
    )
    return selection_group


def parse_input(parser):
    """Parse the input file argument."""
    parser.add_argument(
        "--input",
        "-i",
        type=filepath,
        required=True,
        help="Input file path",
    )


def parse_output(parser, required=True):
    """Parse the output file argument."""
    parser.add_argument(
        "--output",
        "-o",
        type=filepath,
        required=required,
        help="Output file path",
    )


def parse_test(parser):
    """Parse the test flag argument"""
    parser.add_argument("--test", action="store_true", help="Run in test mode")


def create_tasks_parser(subparsers):
    """Create the argument parser for ``openproblems-cli tasks``."""
    subparsers.add_parser("tasks", help="List tasks")


def create_list_parser(subparsers):
    """Create the argument parser for ``openproblems-cli list``."""
    parser = subparsers.add_parser("list", help="List datasets, methods and metrics")
    parser.add_argument(
        "--task",
        "-t",
        type=str,
        help="List methods from a specific task",
    )
    parse_function_type(parser)


def create_image_parser(subparsers):
    """Create the argument parser for ``openproblems-cli image``."""
    parser = subparsers.add_parser(
        "image", help="Fetch a Docker image associated with a function"
    )
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        help="Select functions from a specific task",
        required=True,
    )
    parse_function_type(parser)
    parser.add_argument("name", type=str, help="Name of the selected method")


def create_hash_parser(subparsers):
    """Create the argument parser for ``openproblems-cli hash``."""
    parser = subparsers.add_parser(
        "hash", help="Fetch a git hash associated with a function"
    )
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        help="Select functions from a specific task",
        required=True,
    )
    parse_function_type(parser)
    parser.add_argument("name", type=str, help="Name of the selected method")


def create_load_parser(subparsers):
    """Create the argument parser for ``openproblems-cli load``."""
    parser = subparsers.add_parser("load", help="Load a dataset for a task")
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        help="Select datasets from a specific task",
        required=True,
    )
    parse_test(parser)
    parse_output(parser)
    parser.add_argument("name", type=str, help="Name of the selected dataset")


def create_run_parser(subparsers):
    """Create the argument parser for ``openproblems-cli run``."""
    parser = subparsers.add_parser("run", help="Run a method on a dataset")
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        help="Select methods from a specific task",
        required=True,
    )
    parse_input(parser)
    parse_output(parser)
    parser.add_argument(
        "--version-file",
        type=str,
        help="Write method version to a file",
        default=None,
    )
    parser.add_argument("name", type=str, help="Name of the selected method")
    parse_test(parser)


def create_evaluate_parser(subparsers):
    """Create the argument parser for ``openproblems-cli evaluate``."""
    parser = subparsers.add_parser("evaluate", help="Evaluate a metric on a method")
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        help="Select metrics from a specific task",
        required=True,
    )
    parse_input(parser)
    parse_output(parser, required=False)
    parser.add_argument("name", type=str, help="Name of the selected metric")


def create_test_parser(subparsers):
    """Create the argument parser for ``openproblems-cli test``."""
    parser = subparsers.add_parser("test", help="Test a dataset, method and/or metric")
    parser.add_argument(
        "--task",
        "-t",
        type=str,
        help="Name of task to test",
        required=True,
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Name of the dataset to test. If missing, uses the sample dataset",
        default=None,
    )
    parser.add_argument(
        "--method",
        type=str,
        help="Name of the method to test. If missing, uses the sample method",
        default=None,
    )
    parser.add_argument(
        "--metric",
        type=str,
        help="Name of the metric to test. If missing, no metric is run",
        default=None,
    )
    parse_test(parser)


def create_parser():
    """Create the argument parser for ``openproblems-cli``."""
    parser = argparse.ArgumentParser(
        prog="openproblems-cli",
        description="Open Problems for Single Cell Analysis command-line interface",
    )
    parser.add_argument(
        "--parallel",
        "-p",
        action="store_true",
        help="Run tasks in parallel. This prevents deletion of the cache",
    )
    parser.add_argument(
        "--version",
        "-v",
        action="store_true",
        help="Print version and exit",
    )
    parser.add_argument(
        "--test-hash",
        action="store_true",
        help="Test that `openproblems-cli hash` works",
    )
    subparsers = parser.add_subparsers(dest="subcommand")

    for create_subparser in [
        create_tasks_parser,
        create_list_parser,
        create_image_parser,
        create_load_parser,
        create_run_parser,
        create_evaluate_parser,
        create_hash_parser,
        create_test_parser,
    ]:
        create_subparser(subparsers)

    return parser
