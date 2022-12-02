import subprocess
import os
import json

## VIASH START

par = {
    'input': '/home/kai/Documents/openroblems/openproblems-v2/src/denoising',
    'output': 'output/output.json',
    'show_history': True
}
meta = {
    'functionality_name': 'dca',
}

## VIASH STOP

print(par['show_history'])

output = []

def get_git_file_info(fp, format="none", history=0):
    
    if history:
        cmd = [
        "git",
        'log',
        "--no-merges",
        f"--pretty=format:%H",
        "--",
        fp
    ]
    else:    
        cmd = [
            "git",
            'log',
            "-n",
            "1",
            f"--pretty=format:{format}",
            "--",
            fp
        ]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    return out
    

for root, dirs, files in os.walk(par['input']):
        for file in files:
            git_file= {}
            fp = os.path.join(root,file)
            abs_fp = fp.replace(par['input']+"/","")
            git_file['path'] = abs_fp
            git_file['last_modified'] = get_git_file_info(fp, "%ci")
            git_file['sha'] = get_git_file_info(fp, "%H")
            if "show_history" in par:
                if par['show_history']:
                    git_file['history_sha'] = get_git_file_info(fp,history=1).split("\n")

            output.append(git_file)

with open(par['output'], 'w') as f:
    json.dump(output, f, indent=4)



