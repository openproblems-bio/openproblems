import yaml
from typing import Dict


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

def check_metric(metric: Dict[str, str])  -> str:
    assert "metric_id" in metric, "metric_id not a field"
    assert "metric_name" in metric, f"metric_name not a field in metric['metric_id']"
    assert "min" in metric, f"min not a field in metric['metric_id']"
    assert "max" in metric, f"max not a field in metric['metric_id']"
    assert "maximize" in metric, f"maximize not a field in metric['metric_id']"
    assert isinstance(metric['metric_id'], str), "not a string"
    assert isinstance(metric['metric_name'], str), "not a string"
    assert isinstance(metric['min'], int), "not an int"
    assert isinstance(metric['max'], (int, str)), "not an int or string (+inf)"
    assert isinstance(metric['maximize'], bool), "not a bool"


print("Load config data", flush=True)
with open(meta["config"], "r") as file:
                config = yaml.safe_load(file)

info = config['functionality']['info']

print("check general fields", flush=True)
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"


print("Check info fields", flush=True)
assert "metrics" in info, "metrics not an info field"
for metric in info["metrics"]:
    check_metric(metric)

print("All checks succeeded!", flush=True)
