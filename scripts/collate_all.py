import json
import os
import sys
import openproblems


def main(task_name, dataset_name, input_file, input_dir, output_file):

    task = eval("openproblems.tasks.{}".format(task_name))
    dataset = getattr(task.datasets, dataset_name)

    result = []
    for method in task.METHODS:
        result_file = os.path.join(input_dir, "{}.json".format(method.__name__))
        with open(result_file, "r") as handle:
            result.append(json.load(handle))

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

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(sys.argv[1:])
