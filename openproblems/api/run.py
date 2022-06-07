from ..utils import temporary
from . import utils

import anndata


@temporary(version="1.0")
def _backport_code_version(adata, fun):
    if "method_code_version" not in adata.uns:
        adata.uns["method_code_version"] = fun.metadata["code_version"]
    return adata


def run_method(adata, task_name, function_name, test):
    """Load a dataset for a task."""
    fun = utils.get_function(task_name, "methods", function_name)
    adata = fun(adata, test=test)
    adata = _backport_code_version(adata, fun)
    return adata


def main(args):
    """Run the ``run`` subcommand."""
    adata = anndata.read_h5ad(args.input)
    adata = run_method(adata, args.task, args.name, args.test)
    adata.write_h5ad(args.output)
    return adata.uns["method_code_version"]
