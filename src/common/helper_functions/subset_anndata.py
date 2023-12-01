"""Helper functions related to subsetting AnnData objects based on the file format 
specifications in the .config.vsh.yaml and slot mapping overrides."""

def read_config_slots_info(config_file, slot_mapping = {}):
    """Read the .config.vsh.yaml to find out which output slots need to be copied to which output file.
    
    Arguments:
    config_file -- Path to the .config.vsh.yaml file (required).
    slot_mapping -- Which slots to retain. Must be a dictionary whose keys are the names 
      of the AnnData structs, and values is another dictionary with destination value
      names as keys and source value names as values.
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
        ```
    """
    import yaml
    import re

    # read output spec from yaml
    with open(config_file, "r") as object_name:
        config = yaml.safe_load(object_name)

    output_struct_slots = {}

    # fetch info on which slots should be copied to which file
    for arg in config["functionality"]["arguments"]:
        # argument is an output file with a slot specification
        if arg["direction"] == "output" and arg.get("info", {}).get("slots"):
            object_name = re.sub("--", "", arg["name"])
            
            struct_slots = arg['info']['slots']
            out = {}
            for (struct, slots) in struct_slots.items():
                out_struct = {}
                for slot in slots:
                    # if slot_mapping[struct][slot['name']] exists, use that as the source slot name
                    # otherwise use slot['name']
                    source_slot = slot_mapping.get(struct, {}).get(slot["name"], slot["name"])
                    out_struct[slot["name"]] = source_slot
                out[struct] = out_struct

            output_struct_slots[object_name] = out

    return output_struct_slots

# create new anndata objects according to api spec
def subset_anndata(adata, slot_info):
    """Create new anndata object according to slot info specifications.
    
    Arguments:
    adata -- An AnnData object to subset (required)
    slot_info -- Which slots to retain, typically one of the items in the output of read_config_slots_info.
      Must be a dictionary whose keys are the names of the AnnData structs, and values is another 
      dictionary with destination value names as keys and source value names as values. 
      """
    import pandas as pd
    import anndata as ad

    structs = ["layers", "obs", "var", "uns", "obsp", "obsm", "varp", "varm"]
    kwargs = {}

    for struct in structs:
        slot_mapping = slot_info.get(struct, {})
        data = {dest : getattr(adata, struct)[src] for (dest, src) in slot_mapping.items()}
        if len(data) > 0:
            if struct in ['obs', 'var']:
                data = pd.concat(data, axis=1)
            kwargs[struct] = data
        elif struct in ['obs', 'var']:
            # if no columns need to be copied, we still need an 'obs' and a 'var' 
            # to help determine the shape of the adata
            kwargs[struct] = getattr(adata, struct).iloc[:,[]]

    return ad.AnnData(**kwargs)