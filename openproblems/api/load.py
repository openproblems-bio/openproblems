from . import utils


def load_dataset(task_name, function_name, test):
    """Load a dataset for a task."""
    fun = utils.get_function(task_name, "datasets", function_name)
    return fun(test=test)


def main(args):
    """Run the ``load`` subcommand."""
    adata = load_dataset(args.task, args.name, args.test)
    utils.write_h5ad(adata, args.output)
