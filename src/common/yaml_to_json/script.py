from os import path
import yaml
import json

## VIASH START
par = {
    "input" : ".",
    "task_id" : "denoising",
    "output": "output/task.json",

}
meta = { "functionality" : "foo" }

## VIASH END

with open(par["input"], "r") as f:
    yaml_file = yaml.safe_load(f)


with open(par["output"], "w") as out:
    json.dump(yaml_file, out, indent=2)