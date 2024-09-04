import sys
import torch
from torch.utils.data import DataLoader

import anndata as ad
import pickle
import numpy as np
from scipy.sparse import csc_matrix

#check gpu available
if (torch.cuda.is_available()):
    device = 'cuda:0' #switch to current device
    print('current device: gpu', flush=True)
else:
    device = 'cpu'
    print('current device: cpu', flush=True)


## VIASH START

par = {
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/normal/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/normal/test_mod1.h5ad',
    'input_model': 'resources_test/predict_modality/neurips2021_bmmc_cite/model.pt',
    'input_transform': 'transformer.pickle'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/novel',
    'functionality_name': '171129'
}
## VIASH END

sys.path.append(meta['resources_dir'])
from helper_functions import ModelRegressionAtac2Gex, ModelRegressionAdt2Gex, ModelRegressionGex2Adt, ModelRegressionGex2Atac, ModalityMatchingDataset

print("Load data", flush=True)

input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])

mod1 = input_test_mod1.uns['modality']
mod2 = input_train_mod2.uns['modality']

n_vars_mod1 = input_train_mod2.uns["model_dim"]["mod1"]
n_vars_mod2 = input_train_mod2.uns["model_dim"]["mod2"]

input_test_mod1.X = input_test_mod1.layers['normalized'].tocsr()

# Remove vars that were removed from training set. Mostlyy only applicable for testing.
if input_train_mod2.uns.get("removed_vars"):
  rem_var = input_train_mod2.uns["removed_vars"]
  input_test_mod1 = input_test_mod1[:, ~input_test_mod1.var_names.isin(rem_var)]

del input_train_mod2


model_fp = par['input_model']

print("Start predict", flush=True)

if mod1 == 'GEX' and mod2 == 'ADT':
  model = ModelRegressionGex2Adt(n_vars_mod1,n_vars_mod2)   
  weight = torch.load(model_fp, map_location='cpu')    
  with open(par['input_transform'], 'rb') as f:
    lsi_transformer_gex = pickle.load(f)
  
  model.load_state_dict(weight)    
  input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)

elif mod1 == 'GEX' and mod2 == 'ATAC':
  model = ModelRegressionGex2Atac(n_vars_mod1,n_vars_mod2)   
  weight = torch.load(model_fp, map_location='cpu')
  with open(par['input_transform'], 'rb') as f:
    lsi_transformer_gex = pickle.load(f)
  
  model.load_state_dict(weight)    
  input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)
    
elif mod1 == 'ATAC' and mod2 == 'GEX':
  model = ModelRegressionAtac2Gex(n_vars_mod1,n_vars_mod2)   
  weight = torch.load(model_fp, map_location='cpu')
  with open(par['input_transform'], 'rb') as f:
    lsi_transformer_gex = pickle.load(f)
      
  model.load_state_dict(weight)    
  input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)

elif mod1 == 'ADT' and mod2 == 'GEX':
    model = ModelRegressionAdt2Gex(n_vars_mod1,n_vars_mod2)   
    weight = torch.load(model_fp, map_location='cpu')

    model.load_state_dict(weight)    
    input_test_mod1_ = input_test_mod1.to_df()
    
dataset_test = ModalityMatchingDataset(input_test_mod1_, None, is_train=False)
dataloader_test = DataLoader(dataset_test, 32, shuffle = False, num_workers = 4)

outputs = []
model.eval()
with torch.no_grad():
    for x in dataloader_test:
        output = model(x.float())
        outputs.append(output.detach().cpu().numpy())

outputs = np.concatenate(outputs)
outputs[outputs<0] = 0
outputs = csc_matrix(outputs)

adata = ad.AnnData(
    layers={"normalized": outputs},
    shape=outputs.shape,
    uns={
        'dataset_id': input_test_mod1.uns['dataset_id'],
        'method_id': meta['functionality_name'],
    },
)
adata.write_h5ad(par['output'], compression = "gzip")


