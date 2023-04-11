import yaml
from typing import Dict


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

def check_metric(metric: Dict[str, str])  -> str:
    assert "name" in metric is not None, "name not a field or is empty"
    assert "pretty_name" in metric is not None, "pretty_name not a field in metric or is empty"
    assert "summary" in metric is not None, "summary not a field in metric or is empty"
    assert "description" in metric is not None, "description not a field in metric or is empty"
    assert "reference" in metric, "reference not a field in metric"
    assert "documentation_url" in metric , "documentation_url not a field in metric"
    assert "repository_url" in metric , "repository_url not an info field"
    assert "min" in metric is not None, f"min not a field in metric or is emtpy"
    assert "max" in metric is not None, f"max not a field in metric or is empty"
    assert "maximize" in metric is not None, f"maximize not a field in metric or is emtpy"
    assert isinstance(metric['min'], int), "not an int"
    assert isinstance(metric['max'], (int, str)), "not an int or string (+inf)"
    assert isinstance(metric['maximize'], bool) or metric["maximize"] not in ["-inf", "+inf"], "not a bool"


print("Load config data", flush=True)
with open(meta["config"], "r") as file:
                config = yaml.safe_load(file)

print("check general fields", flush=True)
assert "name" in config["functionality"] is not None, "Name not a field or is empty"
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"


print("Check info fields", flush=True)
info = config['functionality']['info']
print(info)
assert "type" in info, "type not an info field"
assert info["type"] == "metric" , f"got {info['type']} expected 'metric'"
assert "metrics" in info, "metrics not an info field"
for metric in info["metrics"]:
    check_metric(metric)

print("All checks succeeded!", flush=True)
