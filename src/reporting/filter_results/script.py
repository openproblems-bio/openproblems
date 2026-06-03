## VIASH START
par = {
    "input_dataset_info": "resources_test/openproblems/task_results_v4/processed/dataset_info.json",
    "input_method_info": "resources_test/openproblems/task_results_v4/processed/method_info.json",
    "input_metric_info": "resources_test/openproblems/task_results_v4/processed/metric_info.json",
    "input_results": "resources_test/openproblems/task_results_v4/processed/results.json",
    "output_dataset_info": "resources_test/openproblems/task_results_v4/processed/filtered_dataset_info.json",
    "output_method_info": "resources_test/openproblems/task_results_v4/processed/filtered_method_info.json",
    "output_metric_info": "resources_test/openproblems/task_results_v4/processed/filtered_metric_info.json",
    "output_results": "resources_test/openproblems/task_results_v4/processed/filtered_results.json",
    "datasets_exclude": [
        "cellxgene_census/tabula_sapiens",
        "cellxgene_census/mouse_pancreas_atlas",
    ],
    "datasets_include": None,
    "methods_exclude": None,
    "methods_include": None,
    "metrics_exclude": None,
    "metrics_include": None,
}
meta = {"resources_dir": "target/executable/reporting/filter_results"}
## VIASH END

import json
import subprocess
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional


def validate_filtering_args():
    """Validate that include/exclude arguments are mutually exclusive."""
    if par["datasets_include"] and par["datasets_exclude"]:
        raise ValueError(
            "Cannot specify both --datasets_include and --datasets_exclude"
        )

    if par["methods_include"] and par["methods_exclude"]:
        raise ValueError("Cannot specify both --methods_include and --methods_exclude")

    if par["metrics_include"] and par["metrics_exclude"]:
        raise ValueError("Cannot specify both --metrics_include and --metrics_exclude")


def apply_name_filter(
    data_list: List[Dict[str, Any]],
    include_list: Optional[List[str]] = None,
    exclude_list: Optional[List[str]] = None,
    item_type: str = "item",
) -> List[Dict[str, Any]]:
    """Apply filtering to a list based on name field."""
    if not data_list:
        return data_list

    original_count = len(data_list)
    item_names = [item["name"] for item in data_list]

    if include_list:
        items_to_include = set(item_names) & set(include_list)
        if not items_to_include:
            print(
                f"Warning: None of the specified {item_type}s to include were found in the data",
                file=sys.stderr,
            )
            return []

        missing_items = set(include_list) - set(item_names)
        if missing_items:
            print(
                f"Warning: The following {item_type}s specified in include list were not found: "
                + ", ".join(missing_items),
                file=sys.stderr,
            )

        filtered_data = [item for item in data_list if item["name"] in items_to_include]
        print(f"Included {len(filtered_data)} out of {original_count} {item_type}s")
        return filtered_data

    elif exclude_list:
        items_to_exclude = set(item_names) & set(exclude_list)

        missing_items = set(exclude_list) - set(item_names)
        if missing_items:
            print(
                f"Warning: The following {item_type}s specified in exclude list were not found: "
                + ", ".join(missing_items),
                file=sys.stderr,
            )

        filtered_data = [
            item for item in data_list if item["name"] not in items_to_exclude
        ]
        print(
            f"Excluded {len(items_to_exclude)} {item_type}s, keeping {len(filtered_data)} out of {original_count} {item_type}s"
        )
        return filtered_data

    # No filtering applied
    return data_list


def filter_results_data(
    results_data: List[Dict[str, Any]],
    dataset_names: List[str],
    method_names: List[str],
    metric_names: List[str],
) -> List[Dict[str, Any]]:
    """Filter results based on dataset, method, and metric filters."""
    if not results_data:
        return results_data

    original_count = len(results_data)

    # Filter result entries based on dataset_name, method_name, and metric_names
    filtered_results = []
    for result in results_data:
        dataset_keep = result["dataset_name"] in dataset_names
        method_keep = result["method_name"] in method_names

        # Check whether this result should be kept
        if dataset_keep and method_keep:
            filtered_result = result.copy()

            filtered_metrics = [
                (i, name)
                for i, name in enumerate(result["metric_names"])
                if name in metric_names
            ]

            # store metric names
            filtered_result["metric_names"] = [name for _, name in filtered_metrics]

            # store metric values
            filtered_result["metric_values"] = [
                result["metric_values"][i] for i, _ in filtered_metrics
            ]

            # store metric components
            new_metric_components = []
            for component in result.get("metric_components", []):
                new_component = component.copy()
                new_component["metric_names"] = [
                    name for name in component["metric_names"] if name in metric_names
                ]

                # if metric_names are not empty
                if new_component["metric_names"]:
                    new_metric_components.append(new_component)
            filtered_result["metric_components"] = new_metric_components

            filtered_results.append(filtered_result)

    print(
        f"Filtered results: keeping {len(filtered_results)} out of {original_count} result entries"
    )
    return filtered_results


def validate_json_against_schema(
    json_file: str, schema_file: str, name: str
) -> tuple[bool, str]:
    """Validate a JSON file against its schema using ajv-cli.

    Returns:
        tuple[bool, str]: (is_valid, error_message)
    """
    try:
        cmd = [
            "ajv",
            "validate",
            "--spec",
            "draft2020",
            "-s",
            schema_file,
            "-r",
            str(Path(meta["resources_dir"]) / "schemas" / "results_v4" / "core.json"),
            "-d",
            json_file,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✓ {name} validation passed")
            return True, ""
        else:
            error_msg = ""
            if result.stderr:
                error_msg += f"stderr: {result.stderr.strip()}"
            if result.stdout:
                error_msg += f"\nstdout: {result.stdout.strip()}"
            if not error_msg:
                error_msg = "Unknown validation error"

            return False, error_msg

    except FileNotFoundError:
        return False, "ajv-cli not found. Cannot validate schema"


print("====== Filter results ======")

# Validation
print("\n>>> Validating arguments...")
validate_filtering_args()

# Read input files
print("\n>>> Reading input files...")

print(f'Reading dataset info from "{par["input_dataset_info"]}"...')
with open(par["input_dataset_info"], "r") as f:
    dataset_info = json.load(f)

print(f'Reading method info from "{par["input_method_info"]}"...')
with open(par["input_method_info"], "r") as f:
    method_info = json.load(f)

print(f'Reading metric info from "{par["input_metric_info"]}"...')
with open(par["input_metric_info"], "r") as f:
    metric_info = json.load(f)

print(f'Reading results from "{par["input_results"]}"...')
with open(par["input_results"], "r") as f:
    results = json.load(f)

# Apply filters
print("\n>>> Applying filters...")

print("Filtering datasets...")
filtered_dataset_info = apply_name_filter(
    dataset_info, par["datasets_include"], par["datasets_exclude"], "dataset"
)

print("Filtering methods...")
filtered_method_info = apply_name_filter(
    method_info, par["methods_include"], par["methods_exclude"], "method"
)

print("Filtering metrics...")
filtered_metric_info = apply_name_filter(
    metric_info, par["metrics_include"], par["metrics_exclude"], "metric"
)

# Get names for results filtering
filtered_dataset_names = [item["name"] for item in filtered_dataset_info]
filtered_method_names = [item["name"] for item in filtered_method_info]
filtered_metric_names = [item["name"] for item in filtered_metric_info]

print("Filtering results...")
filtered_results = filter_results_data(
    results, filtered_dataset_names, filtered_method_names, filtered_metric_names
)

# Write and validate output files
print("\n>>> Writing and validating output files...")
results_schemas_dir = Path(meta["resources_dir"]) / "schemas" / "results_v4"

validation_files = [
    {
        "data": filtered_dataset_info,
        "schema": "dataset_info.json",
        "file": par["output_dataset_info"],
        "name": "dataset info",
    },
    {
        "data": filtered_method_info,
        "schema": "method_info.json",
        "file": par["output_method_info"],
        "name": "method info",
    },
    {
        "data": filtered_metric_info,
        "schema": "metric_info.json",
        "file": par["output_metric_info"],
        "name": "metric info",
    },
    {
        "data": filtered_results,
        "schema": "results.json",
        "file": par["output_results"],
        "name": "results",
    },
]

all_valid = True
for validation in validation_files:
    print(f'Writing {validation["name"]} to "{validation["file"]}"...')
    with open(validation["file"], "w") as f:
        json.dump(validation["data"], f, indent=2, ensure_ascii=False)

    print(f'Validating {validation["name"]}...')
    schema_file = str(results_schemas_dir / validation["schema"])
    is_valid, error_msg = validate_json_against_schema(
        validation["file"], schema_file, validation["name"]
    )
    if not is_valid:
        print(f'✗ {validation["name"]} validation failed')
        print(f"Validation error: {error_msg}")
        all_valid = False

if not all_valid:
    raise RuntimeError("One or more output files do not conform to their schemas")

# Summary
print("\n>>> Summary of filtering results:")
print(f"Datasets: {len(filtered_dataset_info)} (from {len(dataset_info)})")
print(f"Methods: {len(filtered_method_info)} (from {len(method_info)})")
print(f"Metrics: {len(filtered_metric_info)} (from {len(metric_info)})")
print(f"Results: {len(filtered_results)} (from {len(results)})")

print("\n>>> Done!")
