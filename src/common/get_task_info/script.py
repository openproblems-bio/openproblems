from os import path
from yaml import load, CSafeLoader
import json

## VIASH START
par = {
    "input" : "src",
    "query" : "denoising",
    "output": "output/task.json",

}
meta = { "functionality" : "foo" }

## VIASH END

task_info_path = path.join(par['input'], par['query'], "docs", "task_description.yaml")



with open(task_info_path, 'r') as f:
    task_info = load(f, Loader=CSafeLoader )

with open(par["output"], 'w') as out:
    json.dump(task_info, out, indent=2)