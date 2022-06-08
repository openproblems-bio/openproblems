## VIASH START
par = {
    'input': '../../../data/test_data_preprocessed.h5ad',
    'output': 'output.knnscran.h5ad'
}
## VIASH END
resources_dir = "../../../../common/tools/"
utils_dir = "../../../"

import sys
sys.path.append(resources_dir)
sys.path.append(utils_dir)
sys.path.append(meta['resources_dir'])
import scanpy as sc
from normalize import log_cpm
from utils import classifier
import sklearn.linear_model

print("Load input data")
adata = sc.read(par['input'])

print("Run classifier")
adata = log_cpm(adata)
max_iter = par['max_iter']
adata = classifier(adata, estimator=sklearn.linear_model.LogisticRegression, max_iter=max_iter)
adata.uns["method_id"] = "logistic_regression_log_cpm"

print("Write data")
adata.write(par['output'], compression="gzip")
