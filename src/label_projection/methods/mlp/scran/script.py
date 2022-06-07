## VIASH START
par = {
    'input': '../../../data/test_data_preprocessed.h5ad',
    'output': 'output.mlpscran.h5ad'
}
## VIASH END
resources_dir = "../../../../common/tools/"
utils_dir = "../../../"

import sys
sys.path.append(resources_dir)
sys.path.append(utils_dir)
sys.path.append(meta['resources_dir'])
import scanpy as sc
from normalize import log_scran_pooling
from utils import classifier
import sklearn.neural_network

print("Load input data")
adata = sc.read(par['input'])

print("Run classifier")
hidden_layer_sizes = tuple(par['hidden_layer_sizes'])
max_iter = par['max_iter']
adata = log_scran_pooling(adata)
adata = classifier(adata, estimator=sklearn.neural_network.MLPClassifier, hidden_layer_sizes=hidden_layer_sizes, max_iter=max_iter)
adata.uns["method_id"] = "mlp_scran"

print("Write data")
adata.write(par['output'], compression="gzip")
