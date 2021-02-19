import openproblems
import sys
import json
import anndata


def main(task_name, metric_name, input_file, output_file):
    """Apply a metric to a single method dataset pair."""
    openproblems.data.no_cleanup()
    task = eval("openproblems.tasks.{}".format(task_name))
    metrics = getattr(task, "metrics")
    metric = getattr(metrics, metric_name)
    adata = anndata.read_h5ad(input_file)
    result = metric(adata)
    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(*sys.argv[1:])
