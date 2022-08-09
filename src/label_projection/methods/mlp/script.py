## VIASH START
par = {
    'input': '../../../../resources_test/label_projection/pancreas/toy_preprocessed_data.h5ad',
    'output': 'output.mlplogcpm.h5ad'
}
## VIASH END
resources_dir = "../../../../common/tools/"
utils_dir = "../../"

import sys
sys.path.append(resources_dir)
sys.path.append(utils_dir)
sys.path.append(meta['resources_dir'])
import scanpy as sc
from utils import classifier
import sklearn.neural_network


print("Load input data")
adata = sc.read(par['input'])

if "normalization_method" not in adata.uns:
    print("Warning: trying to predict for a not normalized data")


print("Run classifier")
hidden_layer_sizes = tuple(par['hidden_layer_sizes'])
max_iter = par['max_iter']
adata = classifier(adata, estimator=sklearn.neural_network.MLPClassifier, hidden_layer_sizes=hidden_layer_sizes, max_iter=max_iter)
adata.uns["method_id"] = meta["functionality_name"]

print("Write data")
adata.write(par['output'], compression="gzip")
