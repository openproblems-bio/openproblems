## VIASH START
par = {
    'input': '../../../data/test_data_preprocessed.h5ad',
    'output': 'output.knnscran.h5ad'
}
meta = {
    'resources_dir'
}
## VIASH END
resources_dir = "../../../../common/tools/"
utils_dir = "../../../"

import sys
sys.path.append(meta['resources_dir'])
import scanpy as sc
from utils import classifier
import sklearn.neighbors

print("Load input data")
adata = sc.read(par['input'])

print("Run classifier")
adata = classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)
adata.uns["method_id"] = meta["functionality_name"]

print("Write data")
adata.write(par['output'], compression="gzip")
