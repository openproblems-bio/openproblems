import json
import os
import sys
import openproblems


def main(task_name, meta_file, input_dir, output_file):
    """Collate results from all metrics of a task applied to one method dataset pair."""
    openproblems.data.no_cleanup()

    task = eval("openproblems.tasks.{}".format(task_name))

    with open(meta_file, "r") as handle:
        result = json.load(handle)

    for metric in task.METRICS:
        result_file = os.path.join(input_dir, "{}.metric.json".format(metric.__name__))
        with open(result_file, "r") as handle:
            result[metric.metadata["metric_name"]] = json.load(handle)

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(*sys.argv[1:])
