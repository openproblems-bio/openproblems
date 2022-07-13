## VIASH START
par = {
    'input': '../../../data/test_data_preprocessed.h5ad',
    'output': 'output.knnscran.h5ad'
}
## VIASH END
utils_dir = "../../"

import sys
sys.path.append(utils_dir)
sys.path.append(meta['resources_dir'])
import scanpy as sc
from utils import classifier
import sklearn.linear_model

print("Load input data")
adata = sc.read(par['input'])

print("Run classifier")
max_iter = par['max_iter']
adata = classifier(adata, estimator=sklearn.linear_model.LogisticRegression, max_iter=max_iter)
adata.uns["method_id"] = meta["functionality_name"]

print("Write data")
adata.write(par['output'], compression="gzip")
