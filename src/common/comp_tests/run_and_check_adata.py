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
def check_slots(adata, slot_metadata):
    """Check whether an AnnData file contains all for the required
    slots in the corresponding .info.slots field.
    """
    for struc_name, slot_items in slot_metadata.items():
        struc_dict = getattr(adata, struc_name)
        
        for slot_item in slot_items:
            if slot_item.get("required", False):
                assert slot_item["name"] in struc_dict,\
                    f"File '{arg['value']}' is missing slot .{struc_name}['{slot_item['name']}']"


# read viash config
with open(meta["config"], "r") as stream:
    config = yaml.safe_load(stream)

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

# construct command
cmd = [ meta["executable"] ]
for arg in arguments:
    if arg["type"] == "file":
        cmd.extend([arg["name"], arg["value"]])


print(">> Checking whether input files exist", flush=True)
for arg in arguments:
    if arg["type"] == "file" and arg["direction"] == "input":
        assert path.exists(arg["value"]), f"Input file '{arg['value']}' does not exist"

print(">> Running script as test", flush=True)
# out = subprocess.run(cmd, check=True, capture_output=True, text=True)
subprocess.run(cmd, check=True)

print(">> Checking whether output file exists", flush=True)
for arg in arguments:
    if arg["type"] == "file" and arg["direction"] == "output":
        assert path.exists(arg["value"]), f"Output file '{arg['value']}' does not exist"

print(">> Reading h5ad files and checking formats", flush=True)
adatas = {}
for arg in arguments:
    if arg["type"] == "file":
        print(f"Reading and checking {arg['clean_name']}", flush=True)
        adata = ad.read_h5ad(arg["value"])
        slots = arg["info"]["slots"]

        print(f"  {adata}")

        check_slots(adata, slots)

        adatas[arg["clean_name"]] = adata

print("All checks succeeded!", flush=True)
