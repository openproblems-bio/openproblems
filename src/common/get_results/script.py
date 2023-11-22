import yaml
import json
from pandas import read_csv
from datetime import timedelta
import re
import numpy as np

## VIASH START

par = {
    'input_scores': 'output/v2/batch_integration/output.run_benchmark.output_scores.yaml',
    'input_execution': 'output/v2/batch_integration/trace.txt',
    'methods_meta': 'output/get_info/method_info.json',
    'metrics_meta': 'output/get_info/metric_info.json',
    'task_id': 'batch_integration',
    'output': 'output/temp/results.json'
}

meta = {
}

## VIASH END

def load_meta (meta_file):
    '''
        Load the meta information
    '''
    with open(meta_file, 'r') as f:
        meta = json.load(f)
    return meta

def get_info (id, type,meta_info):
    '''
        Get the information of the method or metric
    '''
    key = type + "_id"
    for info in meta_info:
        if info.get(key) == id:
            return info

def fix_values_scaled(metric_result):
    for i, value in enumerate(metric_result):
        if np.isnan(value):
            metric_result[i] = 0.0
    return metric_result

def fix_nan_scaled(metrics):
    for metric in metrics:
        if np.isnan(metrics[metric]) or np.isinf(metrics[metric]):
            metrics[metric] = 0

    return metrics

def fix_nan (metrics):
    for metric in metrics:
        if np.isnan(metrics[metric]):
            metrics[metric] = "NaN"
        elif np.isneginf(metrics[metric]):
            metrics[metric] = "-Inf"
        elif np.isinf(metrics[metric]):
            metrics[metric] = "Inf"

    return metrics

def organise_score (scores):
    '''
    combine all the metric values into one dictionary per method, dataset and normalization
    '''
    score_temp = {}
    for score in scores:
        score_id = score["dataset_id"] + "_" + score["method_id"] + "_" + score["normalization_id"]
        
        if score.get("metric_values") is None:
            score["metric_values"] = [None] * len(score["metric_ids"])
        # for i, value in enumerate(score["metric_values"]):
            # if np.isnan(value):
            #     score["metric_values"][i] = None
        comb_metric = zip(score["metric_ids"], score["metric_values"])
        score["metric_values"] = dict(comb_metric)
        score["task_id"] = par["task_id"]
        del score["metric_ids"]
        if score_temp.get(score_id) is None:
            score_temp[score_id] = score
        else:
            score_temp[score_id]["metric_values"].update(score["metric_values"])
    
    return score_temp

def normalize_scores (scores, method_info, metric_info):
    """
        Normalize the scores
    """
    
    metric_names=list(set([metric["metric_id"] for metric in metric_info]))

    baseline_methods = [method["method_id"] for method in method_info if method["is_baseline"] ]
    metric_not_maximize = [metric["metric_id"] for metric in metric_info if not metric["maximize"]]
    per_dataset = {}
    for id, score in scores.items():
        if per_dataset.get(score["dataset_id"]) is None:
            per_dataset[score["dataset_id"]] = []
            
        per_dataset[score["dataset_id"]].append(score)

    for id, dataset_results in per_dataset.items():
        for result in dataset_results:
            result["scaled_scores"] = result["metric_values"].copy()

        for metric_name in metric_names:
            metric_values = []
            baseline_values = []
            for result in dataset_results:
                if metric_name in result["metric_values"]:
                    metric_values.append(result["metric_values"][metric_name])
                else:
                    result["metric_values"][metric_name]= float("nan")
                    metric_values.append(0.0)
                if result["method_id"] in baseline_methods:
                    if metric_name in result["metric_values"]:
                        baseline_values.append(result["metric_values"][metric_name])
                    
            baseline_values = fix_values_scaled(baseline_values)
            baseline_values = np.array(baseline_values)
            baseline_min = np.nanmin(baseline_values)
            baseline_range = np.nanmax(baseline_values) - baseline_min
            metric_values = np.array(metric_values)
            metric_values -= baseline_min
            metric_values /= np.where(baseline_range != 0, baseline_range, 1)
            
            if metric_name in metric_not_maximize:
                metric_values = 1 - metric_values
            for result, score in zip(dataset_results,metric_values):
                result["scaled_scores"][metric_name] = score
            
    return per_dataset

def convert_size (df, col):
    '''
        Convert the size to MB and to float type
    '''
    mask_kb = df[col].str.contains("KB")
    mask_mb = df[col].str.contains("MB")
    mask_gb = df[col].str.contains("GB")
    if mask_kb.any():
        df.loc[mask_kb, col] = df.loc[mask_kb, col].str.replace(" KB", "").astype(float)/1024

    if mask_mb.any():
        df.loc[mask_mb, col] = df.loc[mask_mb, col].str.replace(" MB", "").astype(float)

    
    if mask_gb.any():
        df.loc[mask_gb, col] = df.loc[mask_gb, col].str.replace(" GB", "").astype(float)*1024
    return df

def convert_duration(duration_str):
    '''
        Convert the duration to seconds
    '''
    components = duration_str.split(" ")
    hours = 0
    minutes = 0
    seconds = 0
    for component in components:
            milliseconds = 0
    for component in components:
        if "h" in component:
            hours = int(component[:-1])
        elif "ms" in component:
            milliseconds = float(component[:-2])
        elif "m" in component:
            minutes = int(component[:-1])
        elif "s" in component:
            seconds = float(component[:-1])
    duration = timedelta(hours=hours, minutes=minutes, seconds=seconds).total_seconds()
    return duration

def join_trace (traces, result):
    '''
        Join the Seqera (nextflow) trace with the scores
    '''
    trace_dict = {}
    for trace in traces:
        id = trace["name"]
        dataset_id = None
        method_id = None
        match = re.search(r'\((.*?)\)', id)
        id_split = id.split(":")
        if len(id_split)>4:
            method_id = id_split[4]
        if match:
            group = match.group(1)
            split_group = group.split(".")
            if len(split_group)>1:
                dataset_id = split_group[0]
        if dataset_id is not None and method_id is not None:
            dict_id = method_id + "_" + dataset_id
            trace_dict[dict_id] = {
                    "duration_sec": trace["realtime"],
                    "cpu_pct": trace["%cpu"],
                    "peak_memory_mb": trace["peak_vmem"],
                    "disk_read_mb": trace["rchar"],
                    "disk_write_mb": trace["wchar"]
                }
    for score in result:
        search_id = result[score]["method_id"] + "_" + result[score]["dataset_id"]+ "/" + result[score]["normalization_id"]
        if search_id in trace_dict:
            result[score]["resources"] = trace_dict[search_id]
    return result

print('Loading inputs', flush=True)
with open(par['input_scores'], 'r') as f:
    scores = yaml.safe_load(f)
execution = read_csv(par['input_execution'], sep='\t')
method_info = load_meta(par['methods_meta'])
metric_info = load_meta(par['metrics_meta'])

print('Organising scores', flush=True)
org_scores = organise_score(scores)

print('Cleaning execution trace', flush=True)
execution = convert_size(execution, "rchar")
execution = convert_size(execution, "wchar")
execution = convert_size(execution, "peak_vmem")
execution["%cpu"].replace("%", "", regex=True, inplace=True)
execution["realtime"] = execution["realtime"].apply(convert_duration)

print('Joining traces and scores', flush=True)
traces = execution.to_dict(orient="records")
org_scores = join_trace(traces, org_scores)

print('Normalizing scores', flush=True)
org_scores = normalize_scores(org_scores, method_info, metric_info)
# fix NaN en inf
for dataset in org_scores.values():
    for scores in dataset:
        scores["metric_values"] = fix_nan(scores["metric_values"])
        scores["scaled_scores"] = fix_nan_scaled(scores["scaled_scores"])
        scores["mean_score"] = np.array(list(scores["scaled_scores"].values())).mean()

print('Writing results', flush=True)
result = [org_scores[id] for id in org_scores]

result = list(np.concatenate(result).flat)

with open (par['output'], 'w') as f:
    json.dump(result, f, indent=4)