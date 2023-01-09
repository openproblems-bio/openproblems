"""
Schema:

# results/{task.__name__}/{dataset.__name__}.json
{
    "name": dataset.metadata["dataset_name"],
    "data_url": dataset.metadata["data_url"],
    "data_reference": dataset.metadata["data_reference"],
    "headers": {
        "names": [
            "Rank",
            "Name",
            "Metric1 Raw",
            "Metric2 Raw",
            ...,
            "Mean score Scaled",
            "Metric1 Scaled",
            ...,
            "Memory (GB)",
            "Runtime (min)",
            "CPU (%)",
            "Paper",
            "Year",
            "Library"
        ],
        "fixed": ["Name", "Paper", "Library"]
    },
    "results": [
        {
            "Name": method.metadata["method_name"],
            "Paper": method.metadata["paper_name"],
            "Paper URL": method.metadata["paper_url"],
            "Year": method.metadata["year"],
            "Library": method.metadata["code_url"],
            "Implementation": "https://github.com/.../path/to/method.py",
            "Version": method.metadata["method_version"],
            "Runtime (min)": runtime,
            "CPU (%)": cpu,
            "Memory (GB)": memory,
            "Rank": rank,
            "Metric1 Raw": metric1_raw,
            "Metric2 Raw": metric2_raw, .
            ..,
            "Mean score Scaled": mean_score,
            "Metric1 Scaled": metric1,
            ...
        },
        ...
    ]
}
"""
import collections
import copy
import json
import numpy as np
import numpyencoder
import openproblems.api.utils
import os
import pandas as pd
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


def parse_size_to_gb(size):
    """Convert a file size to an integer in GB.

    Example
    -------
    >>> parse_size_to_gb("1000 MB")
    1
    """
    number, unit = [string.strip() for string in size.split()]
    return int(float(number) * size_units[unit]) / size_units["GB"]


time_units = {"s": 1, "m": 60, "h": 3600, "d": 3600 * 24}


def parse_time_to_min(time):
    """Convert a duration to an integer in minutes.

    Example
    -------
    >>> parse_time_to_min("2m 30s")
    2.5
    """
    if " " in time:
        return sum([parse_time_to_min(t) for t in time.split(" ")])
    time = time.strip()
    for unit, value in time_units.items():
        if time.endswith(unit):
            number = float(time.replace(unit, ""))
            return number * value / time_units["m"]


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


def drop_baselines(task_name, dataset_results):
    """Remove baseline methods from dataset results."""
    dataset_results = copy.copy(dataset_results)
    method_names = list(dataset_results.keys())
    n_removed = 0
    for method_name in method_names:
        method = openproblems.api.utils.get_function(
            task_name,
            "methods",
            method_name,
        )
        if method.metadata["is_baseline"]:
            n_removed += 1
            del dataset_results[method_name]

    print(f"Dropped {n_removed} baseline methods")
    return dataset_results


def drop_nan_metrics(dataset_results):
    n_removed = 0
    metric_names = list(list(dataset_results.values())[0]["metrics"].keys())
    for metric_name in metric_names:
        metric_scores = np.array(
            [
                dataset_results[method_name]["metrics"][metric_name]
                for method_name in dataset_results
            ]
        )
        if np.all(np.isnan(metric_scores)):
            n_removed += 1
            for method_name in dataset_results:
                del dataset_results[method_name]["metrics"][metric_name]
                del dataset_results[method_name]["metrics_raw"][metric_name]
    if n_removed > 0:
        print(f"[WARN] Removed {n_removed} all-NaN metrics")
    return dataset_results


def compute_ranking(dataset_results):
    """Rank all methods on a specific dataset."""
    metric_sums = np.zeros(len(dataset_results))
    metric_names = list(dataset_results.values())[0]["metrics"].keys()
    method_names = list(dataset_results.keys())
    for metric_name in metric_names:
        metric_scores = np.array(
            [
                dataset_results[method_name]["metrics"][metric_name]
                for method_name in method_names
            ]
        )
        metric_scores[np.isnan(metric_scores) | np.isneginf(metric_scores)] = 0
        metric_scores[np.isinf(metric_scores)] = 1
        metric_sums += metric_scores

    final_ranking = {
        method_names[method_idx]: rank + 1
        for rank, method_idx in enumerate(np.argsort(metric_sums)[::-1])
    }
    for method_name, metrics_sum in zip(method_names, metric_sums):
        dataset_results[method_name]["mean_score"] = metrics_sum / len(metric_names)
    return dataset_results, final_ranking


def dataset_results_to_json(task_name, dataset_name, dataset_results_raw):
    """Convert the raw dataset results to pretty JSON for web."""
    print(
        f"Formatting {len(dataset_results_raw)} methods for {task_name}.{dataset_name}"
    )
    dataset = openproblems.api.utils.get_function(task_name, "datasets", dataset_name)
    output = dict(
        name=dataset.metadata["dataset_name"],
        data_url=dataset.metadata["data_url"],
        data_reference="https://openproblems.bio/"
        f"bibliography#{dataset.metadata['data_reference']}",
        headers=dict(
            names=["Rank", "Name", "Mean score"], fixed=["Name", "Paper", "Library"]
        ),
        results=list(),
    )
    dataset_results_raw = normalize_scores(task_name, dataset_results_raw)
    dataset_results = drop_baselines(task_name, dataset_results_raw)
    dataset_results = drop_nan_metrics(dataset_results)
    dataset_results, ranking = compute_ranking(dataset_results)
    metric_names = set()
    for method_name, rank in ranking.items():
        method_results = dataset_results[method_name]
        method = openproblems.api.utils.get_function(
            task_name,
            "methods",
            method_name,
        )
        result = {
            "Name": method.metadata["method_name"],
            "Paper": method.metadata["paper_name"],
            "Paper URL": "https://openproblems.bio/"
            f"bibliography#{method.metadata['paper_reference']}",
            "Year": method.metadata["paper_year"],
            "Library": method.metadata["code_url"],
            "Implementation": "https://github.com/openproblems-bio/openproblems/"
            f"blob/main/{method.__module__.replace('.', '/')}.py",
            "Version": method_results["code_version"],
            "Runtime (min)": parse_time_to_min(method_results["realtime"]),
            "CPU (%)": float(method_results["%cpu"].replace("%", "")),
            "Memory (GB)": parse_size_to_gb(method_results["peak_rss"]),
            "Rank": rank,
            "Mean score": method_results["mean_score"],
        }
        result_metrics = {}
        for metric_type in ["metrics_raw", "metrics"]:
            metric_type_name = "Raw" if metric_type == "metrics_raw" else "Scaled"
            for metric_name, metric_result in method_results[metric_type].items():
                metric = openproblems.api.utils.get_function(
                    task_name, "metrics", metric_name
                )
                if np.isnan(metric_result):
                    metric_result = "NaN"
                elif np.isneginf(metric_result):
                    metric_result = "-Inf"
                elif np.isinf(metric_result):
                    metric_result = "Inf"
                metric_name_fmt = f"{metric.metadata['metric_name']} {metric_type_name}"
                result_metrics[metric_name_fmt] = metric_result
                metric_names.add(metric_name_fmt)
        result.update(sorted(result_metrics.items()))
        output["results"].append(result)
    output["headers"]["names"].extend(sorted(list(metric_names)))
    output["headers"]["names"].extend(
        [
            "Memory (GB)",
            "Runtime (min)",
            "CPU (%)",
            "Paper",
            "Year",
            "Library",
        ]
    )
    return output, dataset_results_raw


def results_to_json(results, outdir):
    """Convert the full results to pretty JSON for web."""
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for task_name, task_results in results.items():
        for dataset_name, dataset_results in task_results.items():
            results_dir = os.path.join(outdir, task_name)
            if not os.path.isdir(results_dir):
                os.mkdir(results_dir)
            filename = os.path.join(results_dir, "{}.json".format(dataset_name))
            filename_raw = os.path.join(results_dir, "{}.raw.json".format(dataset_name))
            dataset_results_json, dataset_results_raw = dataset_results_to_json(
                task_name, dataset_name, dataset_results
            )
            with open(filename_raw, "w") as handle:
                dump_json(
                    dataset_results_raw,
                    handle,
                )
            if workflow_utils.task_is_incomplete(
                openproblems.api.utils.str_to_task(task_name)
            ):
                print("Skipping stub task")
            else:
                with open(filename, "w") as handle:
                    dump_json(
                        dataset_results_json,
                        handle,
                    )


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
