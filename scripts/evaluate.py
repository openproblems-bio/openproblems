import numpy as np
import json
import os
import copy

import openproblems
import openproblems.test.utils

RESULTS_DIR = os.path.join("..", "website", "data", "results")


def evaluate_method(task, adata, method):
    output = openproblems.tools.decorators.profile(method)(adata)
    result = {
        "metrics": dict(),
    }
    for metric in task.METRICS:
        result[metric.metadata["metric_name"]] = float(metric(adata))

    del adata
    result["Name"] = method.metadata["method_name"]
    result["Paper"] = method.metadata["paper_name"]
    result["Paper URL"] = method.metadata["paper_url"]
    result["Year"] = method.metadata["paper_year"]
    result["Code"] = method.metadata["code_url"]
    result["Memory (GB)"] = float(output["memory_mb"]) / 1024
    result["Memory leaked (GB)"] = float(output["memory_leaked_mb"]) / 1024
    result["Runtime (min)"] = float(output["runtime_s"]) / 60
    return result


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        pass


def save_result(result, task, dataset_name):
    result = copy.copy(result)
    for i in range(len(result)):
        del result[i]["Memory leaked (GB)"]
    result = {
        "name": dataset_name,
        "headers": {
            "names": ["Rank"]
            + [metric.metadata["metric_name"] for metric in task.METRICS]
            + ["Memory (GB)", "Runtime (min)", "Name", "Paper", "Code", "Year"],
            "fixed": ["Name", "Paper", "Website"],
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
            "{}.json".format(dataset_name),
        ),
        "w",
    ) as handle:
        json.dump(result, handle)


def evaluate_dataset(task, dataset):
    adata = dataset(test=True)
    result = []
    for method in task.METHODS:
        r = evaluate_method(task, adata.copy(), method)
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
        result[i]["Rank"] = int(rank)

    save_result(result, task, dataset.__name__)
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
        json.dump(results, handle)
    return results


if __name__ == "__main__":
    main()
