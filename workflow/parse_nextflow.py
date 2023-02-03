"""
Schema:

# content/benchmarks/{task.__name__}/data/results.json
[
    {
        "task_id": task.__name__,
        "commit_sha": "abc123",
        "method_id": method.__name__,
        "dataset_id": dataset.__name__,
        "submission_time": "1970-01-01 00:00:00.000",
        "code_version": openproblems.__version__,
        "resources": {
            "duration_sec": 100.0,
            "cpu_pct": 100.0,
            "peak_memory_mb": 1000.0,
            "disk_read_mb": 1000.0,
            "disk_write_mb": 1000.0,
        }
        "metric_values": {
            metric.__name__: 1.0,
            ...
        }
        "scaled_scores": {
            metric.__name__: 1.0,
            ...
        },
        "mean_score": 1.0
    },
    ...
]
"""
import collections
import copy
import json
import numpy as np
import numpyencoder
import openproblems.api.utils
import os
import pandas as pd
import pathlib
import sys
import warnings
import workflow_utils


def dump_json(obj, fp):
    """Dump to JSON in a numpy-safe fashion."""
    json.dump(
        obj,
        fp,
        indent=4,
        sort_keys=False,
        separators=(", ", ": "),
        ensure_ascii=False,
        cls=numpyencoder.NumpyEncoder,
    )


size_units = {"B": 1, "KB": 10**3, "MB": 10**6, "GB": 10**9, "TB": 10**12}


def parse_size_to_mb(size):
    """Convert a file size to an integer in MB.

    Example
    -------
    >>> parse_size_to_gb("1 GB")
    1000
    """
    number, unit = [string.strip() for string in size.split()]
    return int(float(number) * size_units[unit]) / size_units["MB"]


time_units = {"s": 1, "m": 60, "h": 3600, "d": 3600 * 24}


def parse_time_to_sec(time):
    """Convert a duration to an integer in seconds.

    Example
    -------
    >>> parse_time_to_min("2m 30s")
    150
    """
    if " " in time:
        return sum([parse_time_to_sec(t) for t in time.split(" ")])
    time = time.strip()
    for unit, value in time_units.items():
        if time.endswith(unit):
            number = float(time.replace(unit, ""))
            return number * value / time_units["s"]


def read_trace(filename):
    """Read the execution trace from nextflow."""
    df = pd.read_csv(filename, sep="\t", index_col=0)
    df = df.loc[df["name"].str.startswith("run_method")]
    df = df.loc[df["exit"].astype(str) == "0"]
    df["tag"] = (
        df["name"]
        .str.replace(r"^.*\(", "", regex=True)
        .str.replace(r"\)$", "", regex=True)
    )
    df["task"] = df["tag"].str.replace(":.*", "", regex=True)
    df["method-dataset"] = (
        df["tag"]
        .str.replace("^.*?:", "", regex=True)
        .str.replace(":.*$", "", regex=True)
    )
    df["method"] = df["method-dataset"].str.replace(r"\-.*$", "", regex=True)
    df["dataset"] = df["method-dataset"].str.replace(r"^.*\-", "", regex=True)
    for k in ["method-dataset", "native_id", "hash", "exit", "status", "tag", "name"]:
        del df[k]
    return df


def parse_trace_to_dict(df):
    """Parse the trace dataframe and convert to dict."""
    print(f"Parsing {df.shape[0]} trace records")
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    for task_name in df["task"].unique():
        df_task = df.loc[df["task"] == task_name]
        print(f"{task_name}: {df_task.shape[0]} records")
        for dataset_name in df_task["dataset"].unique():
            df_dataset = df_task.loc[df_task["dataset"] == dataset_name]
            print(f"{task_name}.{dataset_name}: {df_task.shape[0]} records")
            for _, row in df_dataset.iterrows():
                method_name = row["method"]
                results[task_name][dataset_name][method_name] = row.to_dict()
                results[task_name][dataset_name][method_name]["metrics"] = dict()
                for k in ["task", "dataset", "method"]:
                    del results[task_name][dataset_name][method_name][k]
    return results


def parse_metric_results(results_path, results):
    """Add metric results to the trace output."""
    missing_traces = []
    metric_filenames = os.listdir(os.path.join(results_path, "results/metrics"))
    print(f"Loading {len(metric_filenames)} metric results")
    for filename in sorted(metric_filenames):
        with open(
            os.path.join(results_path, "results/metrics", filename), "r"
        ) as handle:
            result = float(handle.read().strip())
        task_name, dataset_name, method_name, metric_name = filename.replace(
            ".metric.txt", ""
        ).split(".")
        try:
            results[task_name][dataset_name][method_name]["metrics"][
                metric_name
            ] = result
        except KeyError:
            missing_traces.append(filename)
    if len(missing_traces) > 0:
        print("Missing execution trace for metrics: ")
        for filename in missing_traces:
            print("    {}".format(filename))
        warnings.warn(
            "Missing execution traces for {} metrics.".format(len(missing_traces))
        )
    return results


def parse_method_versions(results_path, results):
    """Add method versions to the trace output."""
    missing_traces = []
    for filename in os.listdir(os.path.join(results_path, "results/method_versions")):
        with open(
            os.path.join(results_path, "results/method_versions", filename), "r"
        ) as handle:
            code_version = handle.read().strip()
        task_name, dataset_name, method_name = filename.replace(
            ".method.txt", ""
        ).split(".")

        try:
            results[task_name][dataset_name][method_name]["code_version"] = code_version
        except KeyError:
            missing_traces.append(filename)
    if len(missing_traces) > 0:
        print("Missing execution trace for method code versions: ")
        for filename in missing_traces:
            print("    {}".format(filename))
        warnings.warn(
            "Missing execution traces for {} methods.".format(len(missing_traces))
        )
    return results


def normalize_scores(task_name, dataset_results):
    """Normalize method scores to [0, 1] based on baseline method scores."""
    for method_name in dataset_results:
        # store original unnormalized results
        dataset_results[method_name]["metrics_raw"] = copy.copy(
            dataset_results[method_name]["metrics"]
        )
    metric_names = list(list(dataset_results.values())[0]["metrics"].keys())

    for metric_name in metric_names:
        try:
            metric = openproblems.api.utils.get_function(
                task_name, "metrics", metric_name
            )
        except openproblems.api.utils.NoSuchFunctionError as e:
            print(f"[WARN] {e}")
            del dataset_results[method_name]["metrics"][metric_name]
            continue
        metric_scores = np.array(
            [
                dataset_results[method_name]["metrics"][metric_name]
                for method_name in dataset_results
            ]
        )
        baseline_methods = []
        for method_name in list(dataset_results.keys()):
            try:
                method = openproblems.api.utils.get_function(
                    task_name,
                    "methods",
                    method_name,
                )
            except openproblems.api.utils.NoSuchFunctionError as e:
                print(f"[WARN] {e}")
                del dataset_results[method_name]
            if method.metadata["is_baseline"]:
                baseline_methods.append(method_name)
        if len(baseline_methods) < 2:
            # just use all methods as a fallback
            baseline_methods = dataset_results.keys()
        baseline_scores = np.array(
            [
                dataset_results[method_name]["metrics"][metric_name]
                for method_name in baseline_methods
            ]
        )
        baseline_min = np.nanmin(baseline_scores)
        baseline_range = np.nanmax(baseline_scores) - baseline_min
        metric_scores -= baseline_min
        metric_scores /= np.where(baseline_range != 0, baseline_range, 1)
        if not metric.metadata["maximize"]:
            metric_scores = 1 - metric_scores
        for method_name, score in zip(dataset_results, metric_scores):
            dataset_results[method_name]["metrics"][metric_name] = score
    return dataset_results


def fix_values(metric_result):
    if np.isnan(metric_result):
        return "NaN"
    if np.isneginf(metric_result):
        return "-Inf"
    if np.isinf(metric_result):
        return "Inf"
    return metric_result


def fix_values_scaled(metric_result):
    if np.isnan(metric_result) or np.isinf(metric_result):
        return 0
    return metric_result


def dataset_results_to_json(task_name, dataset_name, dataset_results):
    dataset_results = normalize_scores(task_name, dataset_results)
    out = []
    for method_name, method_results in dataset_results.items():
        raw = {k: fix_values(v) for k, v in method_results["metrics_raw"].items()}
        scaled = {k: fix_values_scaled(v) for k, v in method_results["metrics"].items()}
        resources = {
            "duration_sec": parse_time_to_sec(method_results["duration"]),
            "cpu_pct": float(method_results["%cpu"].replace("%", "")),
            "peak_memory_mb": parse_size_to_mb(method_results["peak_rss"]),
            "disk_read_mb": parse_size_to_mb(method_results["rchar"]),
            "disk_write_mb": parse_size_to_mb(method_results["wchar"]),
        }
        result = {
            "task_id": task_name,
            "commit_sha": workflow_utils.get_sha(),
            "method_id": method_name,
            "dataset_id": dataset_name,
            "submission_time": method_results["submit"],
            "code_version": method_results["code_version"],
            "resources": resources,
            "metric_values": raw,
            "scaled_scores": scaled,
            "mean_score": np.array(list(scaled.values())).mean(),
        }
        out.append(result)
    return out


def results_to_json(results, outdir: pathlib.Path):
    """Convert the full results to pretty JSON for web."""
    for task_name, task_results in results.items():
        task_results_out = []
        task_dir = outdir.joinpath(task_name, "data")
        task_dir.mkdir(parents=True, exist_ok=True)
        for dataset_name, dataset_results in task_results.items():
            results_dir = os.path.join(outdir, task_name)
            if not os.path.isdir(results_dir):
                os.mkdir(results_dir)
            task_results_out.extend(
                dataset_results_to_json(task_name, dataset_name, dataset_results)
            )
        with open(task_dir.joinpath("results.json"), "w") as handle:
            dump_json(task_results_out, handle)


def main(results_path, outdir):
    """Parse the nextflow output."""
    df = read_trace(
        os.path.join(results_path, "results/pipeline_info/execution_trace.txt")
    )
    results = parse_trace_to_dict(df)
    results = parse_metric_results(results_path, results)
    results = parse_method_versions(results_path, results)
    results_to_json(results, outdir)
    return 0


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
