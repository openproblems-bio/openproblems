import anndata as ad
import xgboost as xgb

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data", flush=True)
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])
input_layer = "normalized"

print("Transform into integers", flush=True)
input_train.obs["label_int"] = input_train.obs["label"].cat.codes
categories = input_train.obs["label"].cat.categories

print("Convert AnnDatas into datasets", flush=True)
xg_train = xgb.DMatrix(input_train.layers[input_layer], label=input_train.obs["label_int"])
xg_test = xgb.DMatrix(input_test.layers[input_layer])

print("Fit on train data", flush=True)
param = {'objective': 'multi:softmax', 'num_class': len(categories)}
watchlist = [(xg_train, "train")]
xgb_op = xgb.train(param, xg_train, evals=watchlist)

print("Predict on test data", flush=True)
pred = xgb_op.predict(xg_test).astype(int)
input_test.obs["label_pred"] = categories[pred]

print("Write output to file", flush=True)
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par['output'], compression="gzip")