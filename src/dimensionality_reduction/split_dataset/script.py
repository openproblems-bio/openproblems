import anndata as ad
import re
import yaml
import random
import pandas as pd
import numpy as np

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output_train': "train.h5ad",
    'output_test': "test.h5ad",
}
meta = {
    "functionality_name": "split_data"
}
## VIASH END

# read the .config.vsh.yaml to find out which output slots need to be copied to which output file
def read_slots(par, meta):
    # read output spec from yaml
    with open(meta["config"], "r") as file:
        config = yaml.safe_load(file)

    output_struct_slots = {}

    # fetch info on which slots should be copied to which file
    for arg in config["functionality"]["arguments"]:
        if re.match("--output_", arg["name"]):
            file = re.sub("--output_", "", arg["name"])
            
            struct_slots = arg['info']['slots']
            out = {}
            for (struct, slots) in struct_slots.items():
                out[struct] = { slot['name'] : slot['name'] for slot in slots }

            output_struct_slots[file] = out

    return output_struct_slots

# create new anndata objects according to api spec
def subset_anndata(adata_sub, slot_info):
    structs = ["layers", "obs", "var", "uns", "obsp", "obsm", "varp", "varm"]
    kwargs = {}

    for struct in structs:
        slot_mapping = slot_info.get(struct, {})
        data = {dest : getattr(adata_sub, struct)[src] for (dest, src) in slot_mapping.items()}
        if len(data) > 0:
            if struct in ['obs', 'var']:
                data = pd.concat(data, axis=1)
            kwargs[struct] = data
        elif struct in ['obs', 'var']:
            # if no columns need to be copied, we still need an 'obs' and a 'var' 
            # to help determine the shape of the adata
            kwargs[struct] = getattr(adata_sub, struct).iloc[:,[]]

    return ad.AnnData(**kwargs)

print(">> Load Data", flush=True)
adata = ad.read_h5ad(par["input"])

print(">> Figuring out which data needs to be copied to which output file", flush=True)
slot_info_per_output = read_slots(par, meta)

print(">> Creating train data", flush=True)
output_train = subset_anndata(
    adata_sub=adata, 
    slot_info=slot_info_per_output['train']
)

print(">> Creating test data", flush=True)
output_test = subset_anndata(
    adata_sub=adata,
    slot_info=slot_info_per_output['test']
)

print(">> Writing", flush=True)
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
