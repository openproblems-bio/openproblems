from . import utils

import anndata


def run_method(adata, task_name, function_name):
    """Load a dataset for a task."""
    fun = utils.get_function(task_name, "methods", function_name)
    return fun(adata)


def main(args):
    """Run the ``run`` subcommand."""
    adata = anndata.read_h5ad(args.input)
    adata = run_method(adata, args.task, args.name)
    adata.write_h5ad(args.output)
