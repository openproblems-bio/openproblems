## VIASH START
meta = {
    "executable": "target/docker/methods/lstm_gru_cnn_ensemble/lstm_gru_cnn_ensemble",
    "config": "target/docker/methods/lstm_gru_cnn_ensemble/.config.vsh.yaml",
    "resources_dir": "resources"
}
## VIASH END

# helper functions
def run_component(cmd):
    import subprocess
    
    print(f">> Running script as test", flush=True)
    out = subprocess.run(cmd)

    assert out.returncode == 0, f"Script exited with an error. Return code: {out.returncode}"

def check_input_files(arguments):
    from os import path
    print(">> Checking whether input files exist", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "input" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Input file '{arg['value']}' does not exist"

def check_output_files(arguments):
    from os import path
    print(">> Checking whether output file exists", flush=True)
    for arg in arguments:
        if arg["type"] == "file" and arg["direction"] == "output" and arg["required"]:
            assert not arg["must_exist"] or path.exists(arg["value"]), f"Output file '{arg['value']}' does not exist"

    print(">> Reading output files and checking formats", flush=True)
    for arg in arguments:
        if arg["type"] != "file" or arg["direction"] != "output":
            continue
        check_format(arg)

def check_format(arg):
    arg_info = arg.get("info") or {}
    if arg["type"] == "file":
        arg_format = arg_info.get("format", {})
        file_type = arg_format.get("type") or arg_info.get("file_type")

        # Tabular data
        if file_type in ["parquet", "csv", "tsv"]:
            import pandas as pd
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            if file_type == "csv":
                df = pd.read_csv(arg["value"])
            elif file_type == "tsv":
                df = pd.read_csv(arg["value"], sep="\t")
            else:
                df = pd.read_parquet(arg["value"])
            print(f"  {df}")

            arg_columns = arg_format.get("columns") or arg_info.get("columns") or []
            check_dataframe(df, arg_columns, f"File '{arg['value']}'") 
        
        # Hierarchical data
        elif file_type == "json":
            import json
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            with open(arg["value"]) as f:
                data = json.load(f)
            print(f"  {type(data).__name__} with {len(data)} entries" if isinstance(data, (dict, list)) else f"  {data}")

            check_dictionary(data, arg)
        elif file_type == "yaml":
            import yaml
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            with open(arg["value"]) as f:
                data = yaml.safe_load(f)
            print(f"  {type(data).__name__} with {len(data)} entries" if isinstance(data, (dict, list)) else f"  {data}")

            check_dictionary(data, arg)
        
        # AnnData/SpatialData
        elif file_type in ["h5ad", "anndata_hdf5"]:
            import anndata as ad
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            adata = ad.read_h5ad(arg["value"])

            print(f"  {adata}")

            check_anndata(adata, arg_format, f"File '{arg['value']}'") 
        elif file_type == "anndata_zarr":
            import anndata as ad
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            store = ad.read_zarr(arg["value"])
            print(f"  {store}")

            check_anndata(store, arg_format, f"File '{arg['value']}'") 
        elif file_type == "spatialdata_zarr":
            import spatialdata
            print(f"Reading and checking {arg['clean_name']}", flush=True)

            sdata = spatialdata.read_zarr(arg["value"])
            print(f"  {sdata}")

            check_spatialdata(sdata, arg)


def check_anndata(adata, format_spec, label=""):
    """Check whether an AnnData object contains all required slots
    defined in the given format spec dict.
    """
    for struc_name, items in format_spec.items():
        # skip metadata fields that are not AnnData slots
        if not hasattr(adata, struc_name):
            continue

        struc_x = getattr(adata, struc_name)
        
        if struc_name == "X":
            if items.get("required", True):
                assert struc_x is not None,\
                    f"{label} is missing slot .{struc_name}"
        
        else:
            for item in items:
                if item.get("required", True):
                    assert item["name"] in struc_x,\
                        f"{label} is missing slot .{struc_name}['{item['name']}']" 

def check_dataframe(df, columns, label=""):
    """Check whether a DataFrame contains all required columns
    defined in the given columns spec list.
    """
    for item in columns:
        if item.get("required", True):
            assert item['name'] in df.columns,\
                f"{label} is missing column '{item['name']}'"

def check_dictionary(data, arg):
    """Check whether a JSON/YAML object contains all required top-level keys
    in the corresponding .info.format.keys field.
    """
    arg_info = arg.get("info") or {}
    arg_format = arg_info.get("format", {})
    arg_keys = arg_format.get("keys") or arg_info.get("keys") or []
    for item in arg_keys:
        if item.get("required", True):
            assert isinstance(data, dict) and item["name"] in data,\
                f"File '{arg['value']}' is missing key '{item['name']}'"

def check_spatialdata(sdata, arg):
    """Check whether a SpatialData object contains all required elements
    in the corresponding .info.format field. Supported element categories:
    images, labels, points, shapes, tables.
    """
    arg_info = arg.get("info") or {}
    arg_format = arg_info.get("format") or {}
    element_categories = ["images", "labels", "points", "shapes", "tables"]
    for category in element_categories:
        items = arg_format.get(category) or []
        category_store = getattr(sdata, category, {})
        for item in items:
            if item.get("required", True):
                assert item["name"] in category_store,\
                    f"File '{arg['value']}' is missing {category}['{item['name']}']"
            
            elem_name = item["name"]
            if elem_name not in category_store:
                continue
            element = category_store[elem_name]

            # Check columns for points and shapes (DataFrame-like)
            if category in ["points", "shapes"]:
                check_dataframe(element, item.get("columns") or [], f"File '{arg['value']}' {category}['{elem_name}']")

            # Check obs/var slots for tables (AnnData-like)
            elif category == "tables":
                check_anndata(element, item, f"File '{arg['value']}' tables['{elem_name}']")

def get_argument_sets(config):
    import re
    
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

if __name__ == "__main__":
    import openproblems
    
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
