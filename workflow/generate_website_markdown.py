import openproblems
import os
import pathlib
import re
import sys

INDEX_TOML_TEMPLATE = """+++
title = "{task_name}"
summary = "{task_summary}"
headless = false
theme = "op"
+++
"""

DATASET_TOML_TEMPLATE = """+++
title = "{dataset_name}"
summary = "{dataset_summary}"
+++
"""

API_PATTERN = re.compile(r"^#.*API$")
HEADING_PATTERN = re.compile(r"^# ")


def write_index_md(task, outdir):
    output_md = INDEX_TOML_TEMPLATE.format(
        task_name=task._task_name, task_summary=task._task_summary
    )
    readme_file = task.__file__.replace("__init__.py", "README.md")
    with open(readme_file, "r") as readme_handle:
        for line in readme_handle:
            if HEADING_PATTERN.match(line):
                # exclude top-level headings
                continue
            if API_PATTERN.match(line):
                # exclude everything after ## API
                break
            output_md += line

    output_file = os.path.join(outdir, "_index.md")
    with open(output_file, "w") as output_handle:
        output_handle.write(output_md)


def write_dataset_md(dataset, outdir):
    output_md = DATASET_TOML_TEMPLATE.format(
        dataset_name=dataset.metadata["dataset_name"],
        dataset_summary=dataset.metadata["dataset_summary"],
    )

    dataset_name = dataset.__name__.split(".")[-1]
    output_file = os.path.join(outdir, f"{dataset_name}.md")
    with open(output_file, "w") as output_handle:
        output_handle.write(output_md)


def main(outdir):
    for task in openproblems.TASKS:
        task_outdir = os.path.join(outdir, task.__name__.split(".")[-1])
        if not os.path.isdir(task_outdir):
            pathlib.Path(task_outdir).mkdir(parents=True, exist_ok=True)
        write_index_md(task, task_outdir)
        for dataset in task.DATASETS:
            write_dataset_md(dataset, task_outdir)


if __name__ == "__main__":
    main(sys.argv[1])
