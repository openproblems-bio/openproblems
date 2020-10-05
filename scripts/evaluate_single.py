import anndata
import json
import os
import sys

import openproblems
import openproblems.test.utils


SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))


def evaluate_metric(adata_file, metric_name, output_file):
    subprocess.call(
        [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "evaluate_metric.py"),
            adata_file,
            metric_name,
            output_file,
        ]
    )
    with open(output_file, "r") as handle:
        result = json.load(handle)
    return result


def evaluate_method(task, adata, method):
    output = openproblems.tools.decorators.profile(method)(adata)
    result = dict()
    with tempfile.TemporaryDirectory() as tempdir:
        adata_file = os.path.join(tempdir, "{}.h5ad".format(metric.__name__))
        adata.write_h5ad(adata_file)
        result = []
        for metric in task.METRICS:
            output_file = os.path.join(tempdir, "result.json")
            result[metric.metadata["metric_name"]] = evaluate_metric(
                adata_file, ".".join([metric.__module__, metric.__name__]), output_file
            )

    del adata
    result["Name"] = method.metadata["method_name"]
    result["Paper"] = method.metadata["paper_name"]
    result["Paper URL"] = method.metadata["paper_url"]
    result["Year"] = method.metadata["paper_year"]
    result["Code"] = method.metadata["code_url"]
    result["Version"] = method.metadata["code_version"]
    result["Memory (GB)"] = float(output["memory_mb"] / 1024)
    result["Memory leaked (GB)"] = float(output["memory_leaked_mb"] / 1024)
    result["Runtime (min)"] = float(output["runtime_s"] / 60)
    return result


def main():
    openproblems.test.utils.ignore_numba_warnings()

    task_name = sys.argv[1]
    method_name = sys.argv[2]
    adata_file = sys.argv[3]
    output_file = sys.argv[4]

    task = eval(task_name)
    method = eval(method_name)
    adata = anndata.read_h5ad(adata_file)

    result = evaluate_method(task, adata, method)

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main()
