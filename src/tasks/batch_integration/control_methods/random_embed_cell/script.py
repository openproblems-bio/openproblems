from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import anndata as ad
import yaml

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
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
input.obsm['X_emb'] = OneHotEncoder().fit_transform(
    LabelEncoder().fit_transform(input.obs["label"])[:, None]
)

print("Store outputs", flush=True)
input.uns['output_type'] = output_type
input.uns['hvg'] = par['hvg']
input.uns['method_id'] = meta['functionality_name']
input.write_h5ad(par['output'], compression='gzip')
