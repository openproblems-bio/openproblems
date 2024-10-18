import json
import numpy as np

## VIASH START
## VIASH END

EXPECTED_TASK_FIELDS = ["task_id", "task_name", "task_summary", "task_description"]
EXPECTED_METHOD_FIELDS = ["task_id", "commit_sha", "method_id", "method_name", "method_summary", "paper_reference", "is_baseline"]
EXPECTED_METRIC_FIELDS = ["task_id", "commit_sha", "metric_id", "metric_name", "metric_summary", "paper_reference", "maximize"]
EXPECTED_DATASET_FIELDS = ["task_id", "dataset_id", "dataset_name", "dataset_summary", "data_reference", "data_url"]

def dump_json(obj, fp):
    """Dump to JSON in a numpy-safe fashion."""
    json.dump(
        obj,
        fp,
        indent=4,
        sort_keys=False,
        separators=(", ", ": "),
        ensure_ascii=False,
    )

def create_quality_control(task_info, dataset_info, method_info, metric_info, results):
    """Quality control to detect anomalies in the results."""
    task_id = task_info["task_id"]

    result_qc = []

    def add_qc(
        category: str,
        name: str,
        value,
        severity_value: float,
        code: str,
        message: str,
    ) -> None:
        "Add an entry to the result qc"
        if severity_value <= 1:
            severity = 0
        elif severity_value <= 2:
            severity = 1
        elif severity_value <= 3:
            severity = 2
        else:
            severity = 3
        result_qc.append({
            "task_id": task_id,
            "category": category,
            "name": name,
            "value": value,
            "severity": severity,
            "severity_value": severity_value,
            "code": code,
            "message": message
        })
    
    def percent_missing(list_of_dicts, field):
        are_missing = []
        for item in list_of_dicts:
            if field == 'paper_reference' and item.get('is_baseline', False):
                are_missing.append(0.0)
            elif field in item and item[field] is not None:
                are_missing.append(0.0)
            else:
                are_missing.append(1.0)
        return np.mean(are_missing)
    
    # check task_info
    for field in EXPECTED_TASK_FIELDS:
        pct_missing = percent_missing([task_info], field)
        add_qc(
            "Task info",
            f"Pct '{field}' missing",
            pct_missing,
            3.0 if pct_missing > 0 else 0.0,
            "percent_missing([task_info], field)",
            f"Task metadata field '{field}' should be defined\n"
            f"  Task id: {task_id}\n"
            f"  Field: {field}\n"
        )
    
    # check method_info
    for field in EXPECTED_METHOD_FIELDS:
        pct_missing = percent_missing(method_info, field)
        add_qc(
            "Method info",
            f"Pct '{field}' missing",
            pct_missing,
            3.0 if pct_missing > 0 else 0.0,
            "percent_missing(method_info, field)",
            f"Method metadata field '{field}' should be defined\n"
            f"  Task id: {task_id}\n"
            f"  Field: {field}\n"
        )

    # check metric_info
    for field in EXPECTED_METRIC_FIELDS:
        pct_missing = percent_missing(metric_info, field)
        add_qc(
            "Metric info",
            f"Pct '{field}' missing",
            pct_missing,
            3.0 if pct_missing > 0 else 0.0,
            "percent_missing(metric_info, field)",
            f"Metric metadata field '{field}' should be defined\n"
            f"  Task id: {task_id}\n"
            f"  Field: {field}\n"
        )

    # check dataset_info
    for field in EXPECTED_DATASET_FIELDS:
        pct_missing = percent_missing(dataset_info, field)
        add_qc(
            "Dataset info",
            f"Pct '{field}' missing",
            pct_missing,
            3.0 if pct_missing > 0 else 0.0,
            "percent_missing(dataset_info, field)",
            f"Dataset metadata field '{field}' should be defined\n"
            f"  Task id: {task_id}\n"
            f"  Field: {field}\n"
        )

    # turn results into long format for easier processing
    results_long = [
        {
            "task_id": x["task_id"],
            "method_id": x["method_id"],
            "dataset_id": x["dataset_id"],
            "metric_id": metric["metric_id"],
            "metric_value" : x["metric_values"].get(metric["metric_id"]),
            "scaled_score" : x["scaled_scores"].get(metric["metric_id"]),
        }
        for metric in metric_info
        for x in results
    ]

    # check percentage missing
    pct_missing = 1 - len(results_long) / (len(method_info) * len(metric_info) * len(dataset_info))
    add_qc(
        "Raw data",
        "Number of results",
        len(results),
        pct_missing / .1,
        "len(results) == len(method_info) * len(metric_info) * len(dataset_info)",
        f"Number of results should be equal to #methods × #metrics × #datasets.\n"
        f"  Task id: {task_id}\n"
        f"  Number of results: {len(results)}\n"
        f"  Number of methods: {len(method_info)}\n"
        f"  Number of metrics: {len(metric_info)}\n"
        f"  Number of datasets: {len(dataset_info)}\n"
    )

    # QC per metric
    for metric in metric_info:
        metric_id = metric["metric_id"]
        values = [
            res
            for res in results_long
            if res["metric_id"] == metric_id
            and res["metric_value"] is not None
            and np.isreal(res["metric_value"])
        ]
        pct_missing = 1 - len(values) / len(dataset_info) / len(method_info)

        add_qc(
            "Raw results",
            f"Metric '{metric_id}' %missing",
            pct_missing,
            pct_missing / .1,
            "pct_missing <= .1",
            f"Percentage of missing results should be less than 10%.\n"
            f"  Task id: {task_id}\n"
            f"  Metric id: {metric_id}\n"
            f"  Percentage missing: {pct_missing*100:.0f}%\n"
        )

    # QC per method
    for method in method_info:
        method_id = method["method_id"]
        values = [ 
            res
            for res in results_long
            if res["method_id"] == method_id
            and res["metric_value"] is not None
            and np.isreal(res["metric_value"])
        ]
        pct_missing = 1 - len(values) / len(dataset_info) / len(metric_info)

        add_qc(
            "Raw results",
            f"Method '{method_id}' %missing",
            pct_missing,
            pct_missing / .1,
            "pct_missing <= .1",
            f"Percentage of missing results should be less than 10%.\n"
            f"  Task id: {task_id}\n"
            f"  method id: {method_id}\n"
            f"  Percentage missing: {pct_missing*100:.0f}%\n"
        )

    # QC per dataset
    for dataset in dataset_info:
        dataset_id = dataset["dataset_id"]
        values = [
            res
            for res in results_long
            if res["dataset_id"] == dataset_id
            and res["metric_value"] is not None
            and np.isreal(res["metric_value"])
        ]
        pct_missing = 1 - len(values) / len(metric_info) / len(method_info)

        add_qc(
            "Raw results",
            f"Dataset '{dataset_id}' %missing",
            pct_missing,
            pct_missing / .1,
            "pct_missing <= .1",
            f"Percentage of missing results should be less than 10%.\n"
            f"  Task id: {task_id}\n"
            f"  dataset id: {dataset_id}\n"
            f"  Percentage missing: {pct_missing*100:.0f}%\n"
        )


    # QC per metric and method
    for metric in metric_info:
        for method in method_info:
            metric_id = metric["metric_id"]
            method_id = method["method_id"]
            scores = [ 
                res["scaled_score"]
                for res in results_long
                if res["metric_id"] == metric_id
                and res["method_id"] == method_id
                and res["scaled_score"] is not None
                and np.isreal(res["scaled_score"])
            ]

            if len(scores) >= 1:
                worst_score = np.min(scores).item()
                best_score = np.max(scores).item()

                add_qc(
                    "Scaling",
                    f"Worst score {method_id} {metric_id}",
                    worst_score,
                    worst_score / -1,
                    "worst_score >= -1",
                    f"Method {method_id} performs much worse than baselines.\n"
                    f"  Task id: {task_id}\n"
                    f"  Method id: {method_id}\n"
                    f"  Metric id: {metric_id}\n"
                    f"  Worst score: {worst_score}%\n"
                )

                add_qc(
                    "Scaling",
                    f"Best score {method_id} {metric_id}",
                    best_score,
                    best_score / 2,
                    "best_score <= 2",
                    f"Method {method_id} performs a lot better than baselines.\n"
                    f"  Task id: {task_id}\n"
                    f"  Method id: {method_id}\n"
                    f"  Metric id: {metric_id}\n"
                    f"  Best score: {best_score}%\n"
                )

    return result_qc

def main(par):
    # read data from files
    with open(par["task_info"], "r", encoding="utf8") as file:
        task_info = json.load(file)
    with open(par["method_info"], "r", encoding="utf8") as file:
        method_info = json.load(file)
    with open(par["metric_info"], "r", encoding="utf8") as file:
        metric_info = json.load(file)
    with open(par["dataset_info"], "r", encoding="utf8") as file:
        dataset_info = json.load(file)
    with open(par["results"], "r", encoding="utf8") as file:
        results = json.load(file)

    # create info objects
    quality_control = create_quality_control(task_info, dataset_info, method_info, metric_info, results)

    # write data to files
    with open(par["output"], "w", encoding="utf8") as file:
        dump_json(quality_control, file)

if __name__ == "__main__":
    main(par)
