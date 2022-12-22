import subprocess
import json

## VIASH START
par = {
    'input': '.',
    'output': 'resources_test/input_git_sha.json',
    'show_history': True
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

#? to do: what to do with untracked files?

output = []

def git_ls_files(directory):
    cmd = ["git", "ls-files"]
    cmd_out = subprocess.run(cmd, capture_output=True, text=True, cwd=directory).stdout
    out = [ line for line in cmd_out.split("\n") if line != "" ]
    return out

def get_git_file_info(file, full_history=False):
    # construct command
    cmd = ["git", "log", "--no-merges", "--pretty=%H\t%ci"]

    if not full_history:
        cmd.extend(["-n", "1"])

    cmd.extend(["--", file])
    
    # run command
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=par["input"]).stdout

    # split output
    split = [line.split("\t") for line in out.split("\n") if "\t" in line]

    return split

for relative_path in git_ls_files(par['input']):

    # get git file info
    git_file_info = get_git_file_info(relative_path, full_history=par["show_history"])
    last = git_file_info[0]
    out = {
        "path": relative_path,
        "last_modified": last[1],
        "sha": last[0]
    }
    if par['show_history']:
        out['history_sha'] = [info[0] for info in git_file_info]

    output.append(out)

with open(par['output'], 'w') as f:
    json.dump(output, f, indent=2)



