import collections
import json
import numpy as np
import numpyencoder
import openproblems.api.utils
import os
import pandas as pd
import sys
import utils
import warnings


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
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    for task_name in df["task"].unique():
        df_task = df.loc[df["task"] == task_name]
        for dataset_name in df_task["dataset"].unique():
            df_dataset = df_task.loc[df_task["dataset"] == dataset_name]
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
    for filename in os.listdir(os.path.join(results_path, "results/metrics")):
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


def compute_ranking(task_name, dataset_results):
    """Rank all methods on a specific dataset."""
    rankings = np.zeros(len(dataset_results))
    metric_names = list(dataset_results.values())[0]["metrics"].keys()
    for metric_name in metric_names:
        metric = openproblems.api.utils.get_function(task_name, "metrics", metric_name)
        sorted_order = np.argsort(
            [
                dataset_results[method_name]["metrics"][metric_name]
                for method_name in dataset_results
            ]
        )
        if metric.metadata["maximize"]:
            sorted_order = sorted_order[::-1]
        rankings += np.argsort(sorted_order)

    method_names = list(dataset_results.keys())
    final_ranking = {
        method_names[method_idx]: rank + 1
        for method_idx, rank in zip(
            np.argsort(rankings), np.arange(len(dataset_results))
        )
    }
    return final_ranking


def dataset_results_to_json(task_name, dataset_name, dataset_results):
    """Convert the raw dataset results to pretty JSON for web."""
    dataset = openproblems.api.utils.get_function(task_name, "datasets", dataset_name)
    output = dict(
        name=dataset.metadata["dataset_name"],
        data_url=dataset.metadata["data_url"],
        data_reference=dataset.metadata["data_reference"],
        headers=dict(names=["Rank"], fixed=["Name", "Paper", "Website", "Code"]),
        results=list(),
    )
    ranking = compute_ranking(task_name, dataset_results)
    metric_names = set()
    for method_name, rank in ranking.items():
        method_results = dataset_results[method_name]
        method = openproblems.api.utils.get_function(task_name, "methods", method_name)
        result = {
            "Name": method.metadata["method_name"],
            "Paper": method.metadata["paper_name"],
            "Paper URL": method.metadata["paper_url"],
            "Year": method.metadata["paper_year"],
            "Library": method.metadata["code_url"],
            "Implementation": "https://github.com/openproblems-bio/openproblems/"
            f"blob/main/{method.__module__.replace('.', '/')}",
            "Version": method_results["code_version"],
            "Runtime (min)": parse_time_to_min(method_results["realtime"]),
            "CPU (%)": float(method_results["%cpu"].replace("%", "")),
            "Memory (GB)": parse_size_to_gb(method_results["peak_rss"]),
            "Rank": rank,
        }
        for metric_name, metric_result in method_results["metrics"].items():
            metric = openproblems.api.utils.get_function(
                task_name, "metrics", metric_name
            )
            result[metric.metadata["metric_name"]] = metric_result
            metric_names.add(metric.metadata["metric_name"])
        output["results"].append(result)
    output["headers"]["names"].extend(list(metric_names))
    output["headers"]["names"].extend(
        [
            "Memory (GB)",
            "Runtime (min)",
            "CPU (%)",
            "Name",
            "Paper",
            "Code",
            "Year",
        ]
    )
    return output


def results_to_json(results, outdir):
    """Convert the full results to pretty JSON for web."""
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for task_name, task_results in results.items():
        if utils.task_is_stub(openproblems.api.utils.str_to_task(task_name)):
            # don't write results for stub tasks
            continue
        for dataset_name, dataset_results in task_results.items():
            results_dir = os.path.join(outdir, task_name)
            if not os.path.isdir(results_dir):
                os.mkdir(results_dir)
            filename = os.path.join(results_dir, "{}.json".format(dataset_name))
            try:
                dataset_results_json = dataset_results_to_json(
                    task_name, dataset_name, dataset_results
                )
            except openproblems.api.utils.NoSuchFunctionError:
                continue
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
