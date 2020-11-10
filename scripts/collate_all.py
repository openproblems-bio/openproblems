import json
import os
import sys
import openproblems


def main(input_dir, output_file):
    """Collate results from all tasks."""
    openproblems.data.no_cleanup()

    result = dict()
    for task in openproblems.TASKS:
        task_name = task.__name__.split(".")[-1]
        result[task_name] = dict()
        for dataset in task.DATASETS:
            result[task_name][dataset.__name__] = list()
            for method in task.METHODS:
                result_file = os.path.join(
                    input_dir,
                    task_name,
                    dataset.__name__,
                    "{}.result.json".format(method.__name__),
                )
                with open(result_file, "r") as handle:
                    result[task_name][dataset.__name__].append(json.load(handle))

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(*sys.argv[1:])
