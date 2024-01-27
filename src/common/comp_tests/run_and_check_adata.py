import anndata as ad
import subprocess
from os import path
import yaml
import re

## VIASH START
meta = {
    "executable": "target/docker/denoising/methods/dca/dca",
    "config": "target/docker/denoising/methods/dca/.config.vsh.yaml",
    "resources_dir": "resources_test/denoising"
}
## VIASH END

# helper functions
def check_slots(adata, arg):
    """Check whether an AnnData file contains all for the required
    slots in the corresponding .info.slots field.
    """
    for struc_name, slot_items in arg["info"].get("slots", {}).items():
        struc_x = getattr(adata, struc_name)
        
        if struc_name == "X":
            if slot_items.get("required", True):
                assert struc_x is not None,\
                    f"File '{arg['value']}' is missing slot .{struc_name}"
        
        else:
            for slot_item in slot_items:
                if slot_item.get("required", True):
                    assert slot_item["name"] in struc_x,\
                        f"File '{arg['value']}' is missing slot .{struc_name}['{slot_item['name']}']"

def run_and_check(arguments, cmd):
    print(">> Checking whether input files exist", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "input":
            assert path.exists(arg["value"]), f"Input file '{arg['value']}' does not exist"

    print(f">> Running script as test", flush=True)
    out = subprocess.run(cmd, stderr=subprocess.STDOUT)

    if out.stdout:
        print(out.stdout)

    if out.returncode:
        print(f"script: \'{' '.join(cmd)}\' exited with an error.")
        exit(out.returncode)

    print(">> Checking whether output file exists", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "output":
            assert path.exists(arg["value"]), f"Output file '{arg['value']}' does not exist"

    print(">> Reading h5ad files and checking formats", flush=True)
    adatas = {}
    for arg in arguments:
        if arg["type"] == "file" and "slots" in arg["info"]:
            print(f"Reading and checking {arg['clean_name']}", flush=True)
            adata = ad.read_h5ad(arg["value"])

            print(f"  {adata}")

            check_slots(adata, arg)

            adatas[arg["clean_name"]] = adata

    print("All checks succeeded!", flush=True)


# read viash config
with open(meta["config"], "r") as file:
    config = yaml.safe_load(file)

# get resources
arguments = []

for arg in config["functionality"]["arguments"]:
    new_arg = arg.copy()

    # set clean name
    clean_name = re.sub("^--", "", arg["name"])
    new_arg["clean_name"] = clean_name

    # use example to find test resource file
    if arg["type"] == "file":
      if arg["direction"] == "input":
          value = f"{meta['resources_dir']}/{arg['example'][0]}"
      else:
          value = f"{clean_name}.h5ad"
      new_arg["value"] = value
    
    arguments.append(new_arg)


if "test_setup" not in config["functionality"]["info"]:
    argument_sets = {"run": arguments}
else:
    test_setup = config["functionality"]["info"]["test_setup"]
    argument_sets = {}
    for name, test_instance in test_setup.items():
        new_arguments = []
        for arg in arguments:
            new_arg = arg.copy()
            if arg["clean_name"] in test_instance:
                val = test_instance[arg["clean_name"]]
                if new_arg["type"] == "file" and new_arg["direction"] == "input":
                    val = f"{meta['resources_dir']}/{val}"
                new_arg["value"] = val
            new_arguments.append(new_arg)
        argument_sets[name] = new_arguments

for argset_name, argset_args in argument_sets.items():
    print(f">> Running test '{argset_name}'", flush=True)
    # construct command
    cmd = [ meta["executable"] ]
    for arg in argset_args:
        if arg["type"] == "file":
            cmd.extend([arg["name"], arg["value"]])

    run_and_check(argset_args, cmd)