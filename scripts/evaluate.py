import pandas as pd
import openproblems


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
        result.append(evaluate_method(task, adata, method))

    result = pd.concat(result)
    return result


def evaluate_task(task):
    result = []
    for dataset in task.DATASETS:
        result.append(evaluate_dataset(task, dataset))

    result = pd.concat(result)
    return result


def main():
    result = []
    for task in openproblems.TASKS:
        result.append(evaluate_task(task))

    result = pd.concat(result)
    with open("../results.md", "w") as handle:
        result.to_markdown(handle)
    return result


if __name__ == "__main__":
    main()
