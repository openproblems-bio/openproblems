from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import numpy as np
import anndata as ad
import yaml
from scipy.sparse import csr_matrix

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'jitter': 0.01
}

meta = { 
    'functionality': 'foo',
    'config': 'bar'
}

## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["subtype"]

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])


print('processing data', flush=True)
embedding = OneHotEncoder().fit_transform(
    LabelEncoder().fit_transform(input.obs["label"])[:, None]
)

input.obsm['X_emb'] = csr_matrix(embedding + np.random.uniform(-1 * par['jitter'], par['jitter'], embedding.shape))

print("Store outputs", flush=True)
input.uns['output_type'] = output_type
input.uns['hvg'] = par['hvg']
input.uns['method_id'] = meta['functionality_name']
input.write_h5ad(par['output'], compression='gzip')
