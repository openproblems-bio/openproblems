"""Helper functions related to subsetting AnnData objects based on the file format 
specifications in the .config.vsh.yaml and slot mapping overrides."""

# create new anndata objects according to api spec
def subset_h5ad_by_format(adata, config, arg_name, field_rename_dict = {}):
    """Create new anndata object according to slot info specifications.
    
    Arguments:
    adata -- An AnnData object to subset (required)
    config -- A Viash config object as read by openproblems.project.read_viash_config (required)
    arg_name -- The name of the argument in the config file that specifies the output format (required)
    field_rename_dict -- A mapping between the slots of the source h5ad and the slots of the destination h5ad.
      Example of slot_mapping: 
        ```
        slot_mapping = {
          "layers": {
            "counts": par["layer_counts"],
          },
          "obs": {
            "cell_type": par["obs_cell_type"],
            "batch": par["obs_batch"],
          }
        }
      """
    import pandas as pd
    import anndata as ad

    assert isinstance(adata, ad.AnnData), "adata must be an AnnData object"
    assert isinstance(config, dict), "config must be a dictionary"

    # find argument
    arg = next((x for x in config["all_arguments"] if x["clean_name"] == arg_name), None)
    assert arg, f"Argument '{arg_name}' not found in config"

    # find file format
    file_format = (arg.get("info") or {}).get("format")
    assert file_format, f"Argument '{arg_name}' has no .info.format"

    # find file format type
    file_format_type = file_format.get("type")
    assert file_format_type == "h5ad", "format must be a h5ad type"

    structs = ["layers", "obs", "var", "uns", "obsp", "obsm", "varp", "varm"]
    kwargs = {}

    for struct in structs:
        struct_format = file_format.get(struct, {})
        struct_rename = field_rename_dict.get(struct, {})

        # fetch data from adata
        data = {}
        for field_format in struct_format:
            dest_name = field_format["name"]
            # where to find the data. if the dest_name is in the rename dict, use the renamed name
            # as the source name, otherwise use the dest_name as the source name
            src_name = struct_rename.get(dest_name, dest_name)
            data[dest_name] = getattr(adata, struct)[src_name]
            
        if len(data) > 0:
            if struct in ['obs', 'var']:
                data = pd.concat(data, axis=1)
            kwargs[struct] = data
        elif struct in ['obs', 'var']:
            # if no columns need to be copied, we still need an 'obs' and a 'var' 
            # to help determine the shape of the adata
            kwargs[struct] = getattr(adata, struct).iloc[:,[]]

    return ad.AnnData(**kwargs)
