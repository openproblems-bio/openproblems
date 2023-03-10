import anndata as ad
import scib
import yaml
import re
import pandas as pd

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'hvgs': 2000,
    'output': 'output.h5ad'
}
meta = {}
## VIASH END

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])

def compute_batched_hvg(adata, n_hvgs):
    adata = adata.copy()
    adata.X = adata.layers['normalized'].copy()
    if n_hvgs > adata.n_vars or n_hvgs <= 0:
        hvg_list = adata.var_names.tolist()
    else:
        hvg_list = scib.pp.hvg_batch(
            adata,
            batch_key='batch',
            target_genes=n_hvgs,
            adataOut=False
        )
    adata.var['hvg'] = adata.var_names.isin(hvg_list)
    del adata.X
    return adata


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
            if struct in ["obs", "var"]:
                data = pd.concat(data, axis=1)
            kwargs[struct] = data
        elif struct in ["obs", "var"]:
            # if no columns need to be copied, we still need an "obs" and a "var" 
            # to help determine the shape of the adata
            kwargs[struct] = getattr(adata_sub, struct).iloc[:,[]]

    return ad.AnnData(**kwargs)

print(f'Select {par["hvgs"]} highly variable genes', flush=True)
adata_with_hvg = compute_batched_hvg(input, n_hvgs=par['hvgs'])

print(">> Figuring out which data needs to be copied to which output file", flush=True)
slot_info_per_output = read_slots(par, meta)

print(">> Create output object", flush=True)
output = subset_anndata(
    adata_sub=adata_with_hvg, 
    slot_info=slot_info_per_output["output"]
)

print('Writing adatas to file', flush=True)
output.write(par['output'], compression='gzip')
