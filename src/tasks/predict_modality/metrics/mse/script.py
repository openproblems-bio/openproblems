import anndata as ad
import logging
import numpy as np

## VIASH START
par = {
  "input_test_mod2" : "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  "input_prediction" : "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
  "output" : "openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad"
}
## VIASH END

logging.info("Reading solution file")
ad_sol = ad.read_h5ad(par["input_test_mod2"])

logging.info("Reading prediction file")
ad_pred = ad.read_h5ad(par["input_prediction"])

logging.info("Check prediction format")
if ad_sol.uns["dataset_id"] != ad_pred.uns["dataset_id"]:
  raise ValueError("Prediction and solution have differing dataset_ids")

if ad_sol.shape != ad_pred.shape:
  raise ValueError("Dataset and prediction anndata objects should have the same shape / dimensions.")

logging.info("Computing MSE metrics")

tmp = ad_sol.layers["normalized"] - ad_pred.layers["normalized"]
rmse = np.sqrt(tmp.power(2).mean())
mae = np.abs(tmp).mean()

logging.info("Create output object")
out = ad.AnnData(
  uns = {
    "dataset_id" : ad_pred.uns["dataset_id"],
    "method_id" : ad_pred.uns["method_id"],
    "metric_ids" : ["rmse", "mae"],
    "metric_values" : [rmse, mae],
  }
)

logging.info("Write output to h5ad file")
out.write_h5ad(par["output"], compression=9)
