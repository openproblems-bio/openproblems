import numpy as np
import json
import os
import sys
import tempfile
import subprocess

import openproblems
import openproblems.test.utils

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPTS_DIR, "..", "website", "data", "results")


def evaluate_method(task_name, adata_file, method_name, output_file):
    subprocess.call(
        [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "evaluate_single.py"),
            task_name,
            method_name,
            adata_file,
            output_file,
        ]
    )
    with open(output_file, "r") as handle:
        result = json.load(handle)
    return result


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        pass


def save_result(result, task, dataset):
    result = copy.deepcopy(result)
    for i in range(len(result)):
        del result[i]["Memory leaked (GB)"]
    result = {
        "name": dataset.metadata["dataset_name"],
        "headers": {
            "names": ["Rank"]
            + [metric.metadata["metric_name"] for metric in task.METRICS]
            + ["Memory (GB)", "Runtime (min)", "Name", "Paper", "Code", "Year"],
            "fixed": ["Name", "Paper", "Website", "Code"],
        },
        "results": result,
    }
    results_dir = os.path.join(
        "..", "website", "data", "results", task.__name__.split(".")[-1]
    )
    mkdir(results_dir)
    with open(
        os.path.join(
            results_dir,
            "{}.json".format(dataset.__name__),
        ),
        "w",
    ) as handle:
        json.dump(result, handle, indent=4)


def evaluate_dataset(task, dataset):
    with tempfile.TemporaryDirectory() as tempdir:
        adata = dataset(test=False)
        adata_file = os.path.join(tempdir, "{}.h5ad".format(dataset.__name__))
        adata.write_h5ad(adata_file)
        result = []
        for method in task.METHODS:
            output_file = os.path.join(tempdir, "result.json")
            r = evaluate_method(
                task.__name__,
                adata_file,
                ".".join(method.__module__, method.__name__),
                output_file,
            )
            result.append(r)

    del adata

    rankings = np.zeros(len(result))
    for metric in task.METRICS:
        sorted_order = np.argsort([r[metric.metadata["metric_name"]] for r in result])
        if metric.metadata["maximize"]:
            sorted_order = sorted_order[::-1]
        rankings += np.argsort(sorted_order)

    final_ranking = np.argsort(np.argsort(rankings))
    for i, rank in enumerate(final_ranking):
        result[i]["Rank"] = int(rank) + 1

    save_result(result, task, dataset)
    return result


def evaluate_task(task):
    result = dict()
    for dataset in task.DATASETS:
        result[dataset.__name__] = evaluate_dataset(task, dataset)

    return result


def main():
    openproblems.test.utils.ignore_numba_warnings()

    results = dict()
    mkdir(RESULTS_DIR)
    for task in openproblems.TASKS:
        task_name = task.__name__.split(".")[-1]
        result = evaluate_task(task)
        results[task_name] = result

    with open("../results.json", "w") as handle:
        json.dump(results, handle, indent=4)
    return results


if __name__ == "__main__":
    main()
