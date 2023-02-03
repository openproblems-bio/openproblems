import functools
import git
import json
import openproblems
import pathlib
import re
import sys
import workflow_utils

API_PATTERN = re.compile(r"^#.*API$")
HEADING_PATTERN = re.compile(r"^# ")


@functools.lru_cache()
def get_sha():
    repo = git.Repo(pathlib.Path(openproblems.__path__[0]).parent)
    assert not repo.bare
    return repo.head.commit.hexsha


def get_task_description(task):
    description = ""
    readme_file = task.__file__.replace("__init__.py", "README.md")
    with open(readme_file, "r") as readme_handle:
        for line in readme_handle:
            if HEADING_PATTERN.match(line):
                # exclude top-level headings
                continue
            if API_PATTERN.match(line):
                # exclude everything after ## API
                break
            description += line
    return description


def write_task_json(task, outdir: pathlib.Path):
    data = {
        "id": openproblems.utils.get_member_id(task),
        "commit_sha": get_sha(),
        "task_name": task._task_name,
        "task_summary": task._task_summary,
        "task_description": get_task_description(task),
        "repo": "openproblems-bio/openproblems",
    }
    with open(outdir.joinpath("task_info.json"), "w") as handle:
        json.dump(data, handle, indent=4)


def _write_function_json(task, outdir: pathlib.Path, functions, function_type: str):
    data = []
    for function in functions:
        function.metadata.update(
            {
                "task_id": openproblems.utils.get_member_id(task),
                "commit_sha": get_sha(),
                f"{function_type}_id": openproblems.utils.get_member_id(function),
            }
        )
        data.append(function.metadata)

    with open(outdir.joinpath(f"{function_type}_info.json"), "w") as handle:
        json.dump(data, handle, indent=4)


def write_dataset_json(task, outdir: pathlib.Path):
    _write_function_json(task, outdir, task.DATASETS, "dataset")


def write_method_json(task, outdir: pathlib.Path):
    _write_function_json(task, outdir, task.METHODS, "method")


def write_metric_json(task, outdir: pathlib.Path):
    _write_function_json(task, outdir, task.METRICS, "metric")


def main(outdir: pathlib.Path):
    for task in openproblems.TASKS:
        if workflow_utils.task_is_incomplete(task):
            # don't write json for incomplete tasks
            continue
        task_outdir = outdir.joinpath(openproblems.utils.get_member_id(task), "data")
        task_outdir.mkdir(parents=True, exist_ok=True)
        write_task_json(task, task_outdir)
        write_dataset_json(task, task_outdir)
        write_method_json(task, task_outdir)
        write_metric_json(task, task_outdir)


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
