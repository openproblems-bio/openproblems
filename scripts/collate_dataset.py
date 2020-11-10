import numpy as np
import json
import os
import sys
import openproblems


def main(task_name, dataset_name, input_dir, output_file):
    """Collate results from all methods of a task applied to one dataset."""
    openproblems.data.no_cleanup()

    task = eval("openproblems.tasks.{}".format(task_name))
    dataset = getattr(task.datasets, dataset_name)

    result = []
    for method in task.METHODS:
        result_file = os.path.join(input_dir, "{}.result.json".format(method.__name__))
        with open(result_file, "r") as handle:
            result.append(json.load(handle))

    rankings = np.zeros(len(result))
    for metric in task.METRICS:
        sorted_order = np.argsort([r[metric.metadata["metric_name"]] for r in result])
        if metric.metadata["maximize"]:
            sorted_order = sorted_order[::-1]
        rankings += np.argsort(sorted_order)

    final_ranking = np.argsort(np.argsort(rankings))
    for i, rank in enumerate(final_ranking):
        result[i]["Rank"] = int(rank) + 1
        del result[i]["Memory leaked (GB)"]

    names = ["Rank"]
    names += [metric.metadata["metric_name"] for metric in task.METRICS]
    names += ["Memory (GB)", "Runtime (min)", "Name", "Paper", "Code", "Year"]
    result = {
        "name": dataset.metadata["dataset_name"],
        "headers": {
            "names": names,
            "fixed": ["Name", "Paper", "Website", "Code"],
        },
        "results": result,
    }

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(*sys.argv[1:])
