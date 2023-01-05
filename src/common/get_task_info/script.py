from os import path
from yaml import load, CSafeLoader
import json

## VIASH START
par = {
    "input" : "src",
    "task_id" : "denoising",
    "output": "output/task.json",

}
meta = { "functionality" : "foo" }

## VIASH END

task_info_path = path.join(par['input'], 'src', par['task_id'], "docs", "task_info.yaml")



with open(task_info_path, 'r') as f:
    task_info = load(f, Loader=CSafeLoader )

with open(par["output"], 'w') as out:
    json.dump(task_info, out, indent=2)