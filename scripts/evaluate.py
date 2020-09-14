import json
import os
import openproblems
import openproblems.test.utils

RESULTS_DIR = os.path.join("..", "website", "data", "results")


def evaluate_method(task, adata, method):
    output = openproblems.tools.decorators.profile(method)(adata)
    result = {
        "metrics": dict(),
    }
    for metric in task.METRICS:
        result["metrics"][metric.__name__] = float(metric(adata))

    del adata
    result["method"] = method.metadata["method_name"]
    result["paper_name"] = method.metadata["paper_name"]
    result["paper_url"] = method.metadata["paper_url"]
    result["paper_year"] = method.metadata["paper_year"]
    result["code_url"] = method.metadata["code_url"]
    result["memory_mb"] = float(output["memory_mb"])
    result["memory_leaked_mb"] = float(output["memory_leaked_mb"])
    result["runtime_s"] = float(output["runtime_s"])
    return result


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        pass


def save_result(result, task_name, dataset_name):
    results_dir = os.path.join("..", "website", "data", "results", task_name)
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
        save_result(r, task.__name__.split(".")[-1], dataset.__name__)
        result.append(r)

    del adata
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
