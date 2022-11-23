import re
import yaml
import random
import pandas as pd
import numpy as np
import anndata as ad
# Todo: throw error when not all slots are available?

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'method': 'batch',
    'seed': None,
    'obs_batch': 'batch',
    'obs_label': 'celltype',
    'output_train': 'train.h5ad',
    'output_test': 'test.h5ad',
    'output_solution': 'solution.h5ad'
}
meta = {
    'resources_dir': 'src/label_projection/split_dataset',
    'config': 'src/label_projection/split_dataset/.config.vsh.yaml'
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
            
            # rename source keys
            if 'obs' in out:
                if 'label' in out['obs']:
                    out['obs']['label'] = par['obs_label']
                if 'batch' in out['obs']:
                    out['obs']['batch'] = par['obs_batch']

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

# set seed if need be
if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print(">> Load data")
adata = ad.read_h5ad(par["input"])
print("adata:", adata)

print(f">> Process data using {par['method']} method")
if par["method"] == "batch":
    batch_info = adata.obs[par["obs_batch"]]
    batch_categories = batch_info.dtype.categories
    test_batches = random.sample(list(batch_categories), 1)
    is_test = [ x in test_batches for x in batch_info ]
elif par["method"] == "random":
    train_ix = np.random.choice(adata.n_obs, round(adata.n_obs * 0.8), replace=False)
    is_test = [ not x in train_ix for x in range(0, adata.n_obs) ]

# subset the different adatas
print(">> Figuring which data needs to be copied to which output file")
slot_info_per_output = read_slots(par, meta)

print(">> Creating train data")
output_train = subset_anndata(
    adata_sub=adata[[not x for x in is_test]], 
    slot_info=slot_info_per_output['train']
)

print(">> Creating test data")
output_test = subset_anndata(
    adata_sub=adata[is_test],
    slot_info=slot_info_per_output['test']
)

print(">> Creating solution data")
output_solution = subset_anndata(
    adata_sub=adata[is_test],
    slot_info=slot_info_per_output['solution']
)

print(">> Writing data")
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
output_solution.write_h5ad(par["output_solution"])
