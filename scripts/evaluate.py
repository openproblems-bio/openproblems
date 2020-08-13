import pandas as pd
import openproblems
import openproblems.test.utils


def evaluate_method(task, adata, method):
    method(adata)
    result = {"metric": [], "value": []}
    for metric in task.METRICS:
        result["metric"].append(metric.__name__)
        result["value"].append(metric(adata))

    result = pd.DataFrame(result)
    result["method"] = method.__name__
    result["task"] = task._task_name
    return result


def evaluate_dataset(task, dataset):
    adata = dataset()
    result = []
    for method in task.METHODS:
        r = evaluate_method(task, adata.copy(), method)
        r["dataset"] = dataset.__name__
        result.append(r)

    result = pd.concat(result)
    return result


def evaluate_task(task):
    result = []
    for dataset in task.DATASETS:
        result.append(evaluate_dataset(task, dataset))

    result = pd.concat(result)
    return result


def main():
    openproblems.test.utils.ignore_numba_warnings()

    results = []
    for task in openproblems.TASKS:
        result = evaluate_task(task).sort_values(["dataset", "metric", "value"])
        with open("../website/content/results/{}.md".format(task), "w") as handle:
            result.to_markdown(handle)
        results.append(result)

    results = pd.concat(results).sort_values(["task", "dataset", "metric", "value"])
    with open("../results.md", "w") as handle:
        results.to_markdown(handle)
    return results


if __name__ == "__main__":
    main()
