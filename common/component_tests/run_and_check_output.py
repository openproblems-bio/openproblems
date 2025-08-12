import anndata as ad
import pandas as pd
import subprocess
from os import path
import re
import openproblems

## VIASH START
meta = {
    "executable": "target/docker/methods/lstm_gru_cnn_ensemble/lstm_gru_cnn_ensemble",
    "config": "target/docker/methods/lstm_gru_cnn_ensemble/.config.vsh.yaml",
    "resources_dir": "resources"
}
## VIASH END

# helper functions
def run_component(cmd):
    print(f">> Running script as test", flush=True)
    out = subprocess.run(cmd)

    assert out.returncode == 0, f"Script exited with an error. Return code: {out.returncode}"

def check_input_files(arguments):
    print(">> Checking whether input files exist", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "input" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Input file '{arg['value']}' does not exist"

def check_output_files(arguments):
    print(">> Checking whether output file exists", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "output" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Output file '{arg['value']}' does not exist"

    print(">> Reading h5ad files and checking formats", flush=True)
    for arg in arguments:
        if arg["type"] != "file" or arg["direction"] != "output":
            continue
        check_format(arg)

def check_format(arg):
    arg_info = arg.get("info") or {}
    if arg["type"] == "file":
        arg_format = arg_info.get("format", {})
        file_type = arg_format.get("type") or arg_info.get("file_type")

        if file_type == "h5ad":
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            # try to read as an anndata, else as a parquet file
            adata = ad.read_h5ad(arg["value"])

            print(f"  {adata}")

            check_h5ad_slots(adata, arg)
        elif file_type in ["parquet", "csv", "tsv"]:
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            if file_type == "csv":
                df = pd.read_csv(arg["value"])
            if file_type == "tsv":
                df = pd.read_csv(arg["value"], sep="\t")
            else:
                df = pd.read_parquet(arg["value"])
            print(f"  {df}")
            
            check_df_columns(df, arg)


def check_h5ad_slots(adata, arg):
    """Check whether an AnnData file contains all for the required
    slots in the corresponding .info.format field.
    """
    arg_info = arg.get("info") or {}
    arg_format = arg_info.get("format") or arg_info.get("slots") or {}
    for struc_name, items in arg_format.items():
        # skip the type field
        if struc_name == "type":
            continue

        struc_x = getattr(adata, struc_name)
        
        if struc_name == "X":
            if items.get("required", True):
                assert struc_x is not None,\
                    f"File '{arg['value']}' is missing slot .{struc_name}"
        
        else:
            for item in items:
                if item.get("required", True):
                    assert item["name"] in struc_x,\
                        f"File '{arg['value']}' is missing slot .{struc_name}['{item['name']}']"

def check_df_columns(df, arg):
    """Check whether a DataFrame contains all for the required
    columns in the corresponding .info.columns field.
    """
    arg_info = arg.get("info") or {}
    arg_format = arg_info.get("format", {})
    arg_columns = arg_format.get("columns") or arg_info.get("columns") or []
    for item in arg_columns:
        if item.get("required", True):
            assert item['name'] in df.columns,\
                f"File '{arg['value']}' is missing column '{item['name']}'"

def get_argument_sets(config):
    # get resources
    arguments = []

    for arg in config["all_arguments"]:
        new_arg = arg.copy()
        arg_info = new_arg.get("info") or {}
        example = arg.get("example", [None])[0]

        # use example to find test resource file
        if example and arg["type"] == "file":
            if arg["direction"] == "input":
                value = f"{meta['resources_dir']}/{example}"
            else:
                ext_res = re.search(r"\.(\w+)$", example)
                if ext_res:
                    value = f"{arg['clean_name']}.{ext_res.group(1)}"
                else:
                    value = f"{arg['clean_name']}"
            new_arg["value"] = value
        elif "test_default" in arg_info:
            new_arg["value"] = arg_info["test_default"]
        
        arguments.append(new_arg)

    config_info = config.get("info") or {}
    if "test_setup" not in config_info:
        argument_sets = {"run": arguments}
    else:
        test_setup = config_info["test_setup"]
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

    return argument_sets

def generate_cmd_args(argument_set):
    cmd_args = []
    for arg in argument_set:
        if "value" in arg:
            value = arg["value"]
            if arg["multiple"] and isinstance(value, list):
                value = arg["multiple_sep"].join(value)
            cmd_args.extend([arg["name"], str(value)])
    return cmd_args

# read viash config
config = openproblems.project.read_viash_config(meta["config"])

# get argument sets
argument_sets = get_argument_sets(config)

# run component for each argument set
for argset_name, argset_args in argument_sets.items():
    print(f">> Running test '{argset_name}'", flush=True)
    # construct command
    cmd = [ meta["executable"] ] + generate_cmd_args(argset_args)

    # check input files
    check_input_files(argset_args)

    # run component
    run_component(cmd)

    # check output files
    check_output_files(argset_args)

print("All checks succeeded!", flush=True)