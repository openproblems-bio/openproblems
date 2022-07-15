from . import evaluate
from . import load
from . import run
from . import utils


def _get_api_function(task_name, method_name):
    return utils.get_function(task_name, "api", method_name)


def main(args):
    if args.dataset is not None:
        adata = load.load_dataset(args.task, args.dataset, args.test)
    else:
        adata = _get_api_function(args.task, "sample_dataset")()
    assert _get_api_function(args.task, "check_dataset")(adata)
    if args.method is not None:
        adata = run.run_method(adata, args.task, args.method, args.test)
    else:
        adata = _get_api_function(args.task, "sample_method")(adata)
    assert _get_api_function(args.task, "check_method")(adata)
    if args.metric is not None:
        return evaluate.evaluate_metric(adata, args.task, args.metric)
