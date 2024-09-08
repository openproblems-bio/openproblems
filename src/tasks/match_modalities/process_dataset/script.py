import sys
import random
import numpy as np
import anndata as ad

## VIASH START
par = {
  "input_mod1": "resources_test/common/scicar_cell_lines/dataset_mod1.h5ad",
  "input_mod2": "resources_test/common/scicar_cell_lines/dataset_mod2.h5ad",
  "output_mod1": "output_mod1.h5ad",
  "output_mod2": "output_mod2.h5ad",
  "output_solution_mod1": "output_solution_mod1.h5ad",
  "output_solution_mod2": "output_solution_mod2.h5ad",
  "seed": 123
}
meta = {
  "resources_dir": "src/common/helper_functions/",
  "config": "src/tasks/match_modalities/process_dataset/.config.vsh.yaml"
}
## VIASH END

# import helper functions
sys.path.append(meta["resources_dir"])
from subset_anndata import read_config_slots_info, subset_anndata

# set seed if need be
if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print(">> Load data", flush=True)
input_mod1 = ad.read_h5ad(par["input_mod1"])
input_mod2 = ad.read_h5ad(par["input_mod2"])

print(f">> Permute input data")
mod1_perm = np.random.permutation(np.arange(input_mod1.n_obs))
mod2_perm = np.random.permutation(np.arange(input_mod2.n_obs))

output_mod1 = input_mod1[mod1_perm]
output_mod1.obs_names = [f"cell_mod1_{i}" for i in range(output_mod1.n_obs)]
output_mod2 = input_mod2[mod2_perm]
output_mod2.obs_names = [f"cell_mod2_{i}" for i in range(output_mod2.n_obs)]

print(f">> Create solution objects")
output_solution_mod1 = input_mod1.copy()
output_solution_mod1.obs["permutation_indices"] = np.argsort(mod1_perm)
output_solution_mod2 = input_mod2.copy()
output_solution_mod2.obs["permutation_indices"] = np.argsort(mod2_perm)
    
# subset the different adatas
print(">> Read slot info from config file", flush=True)
slot_info = read_config_slots_info(meta["config"])

print(">> Subset anndatas", flush=True)
output_mod1 = subset_anndata(output_mod1, slot_info["output_mod1"])
output_mod2 = subset_anndata(output_mod2, slot_info["output_mod2"])
output_solution_mod1 = subset_anndata(output_solution_mod1, slot_info["output_solution_mod1"])
output_solution_mod2 = subset_anndata(output_solution_mod2, slot_info["output_solution_mod2"])

print(">> Writing data", flush=True)
output_mod1.write_h5ad(par["output_mod1"])
output_mod2.write_h5ad(par["output_mod2"])
output_solution_mod1.write_h5ad(par["output_solution_mod1"])
output_solution_mod2.write_h5ad(par["output_solution_mod2"])
