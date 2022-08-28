from . import utils

import anndata


def evaluate_metric(adata, task_name, function_name):
    """Load a dataset for a task."""
    fun = utils.get_function(task_name, "metrics", function_name)
    return fun(adata)


def main(args):
    """Run the ``evaluate`` subcommand."""
    adata = anndata.read_h5ad(args.input)
    result = evaluate_metric(adata, args.task, args.name)
    if args.output is not None:
        with open(args.output, "w") as handle:
            handle.write(str(result))
    else:
        return result
