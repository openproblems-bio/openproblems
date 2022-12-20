import json
from typing import Dict, List

## VIASH START
par = {
    'git_sha': 'temp/openproblems-v1.json',
    'comp_info': 'temp/denoising_metrics.yaml',
    'output': 'temp/migration_status.yaml'
}
## VIASH END

def check_status(comp_item: Dict[str, str], git_objects: List[Dict[str, str]]) -> str:
    """Looks for the comp_item's matching git_object 
    based on the comp_item["v1_url"] and git_object["path"].
    If found, checks whether the comp_item["v1_commit"] equals
    git_object["sha"]."""

    v1_url = comp_item.get("v1_url")
    if not v1_url:
        return "v1_url missing"
    
    v1_commit = comp_item.get("v1_commit")
    if not v1_commit:
        return "v1_commit missing"
    
    git_object = [ obj for obj in git_objects if obj["path"] == v1_url ]
    if not git_object:
        return "v1_url does not exist in git repo"

    git_sha = git_object[0]["sha"]
    if git_sha == comp_item["v1_commit"]:
        return "up to date"
    else:
        return f"out of date (sha: {git_sha})"

with open(par['git_sha'], 'r') as f1:
    git_objects = json.load(f1)

with open(par['comp_info'], 'r') as f2:
    comp_items = json.load(f2)

output = []
for comp_item in comp_items:
    # get status
    status = check_status(comp_item, git_objects)

    # store results
    output.append(comp_item | {"status": status})

# write to file
with open(par['output'], 'w') as outf:
    json.dump(output, outf, indent=2)
