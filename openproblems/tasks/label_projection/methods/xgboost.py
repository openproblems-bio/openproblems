from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from typing import Optional

import functools
import numpy as np

_xgboost_method = functools.partial(
    method,
    paper_name="XGBoost: A Scalable Tree Boosting System",
    paper_url="https://doi.org/10.1145/2939672.2939785",
    paper_year=2016,
    code_url="https://xgboost.readthedocs.io/en/stable/index.html",
)


def _xgboost(adata, test: bool = False, num_round: Optional[int] = None):
    import xgboost as xgb

    if test:
        num_round = num_round or 2
    else:  # pragma: nocover
        num_round = num_round or 5

    adata.strings_to_categoricals()
    adata.obs["labels_int"] = adata.obs["labels"].cat.codes
    categories = adata.obs["labels"].cat.categories

    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()

    xg_train = xgb.DMatrix(adata_train.X, label=adata_train.obs["labels_int"])
    xg_test = xgb.DMatrix(adata_test.X, label=adata_test.obs["labels_int"])

    param = dict(
        objective="multi:softmax",
        num_class=len(categories),
    )

    watchlist = [(xg_train, "train")]
    xgb_op = xgb.train(param, xg_train, num_boost_round=num_round, evals=watchlist)

    # Predict on test data
    pred = xgb_op.predict(xg_test).astype(int)
    adata_test.obs["labels_pred"] = categories[pred]

    adata.obs["labels_pred"] = [
        adata_test.obs["labels_pred"][idx] if idx in adata_test.obs_names else np.nan
        for idx in adata.obs_names
    ]

    adata.uns["method_code_version"] = check_version("xgboost")
    return adata


@_xgboost_method(
    method_name="XGBoost (log CPM)",
    image="openproblems-python-extras",
)
def xgboost_log_cpm(adata, test: bool = False, num_round: Optional[int] = None):
    adata = log_cpm(adata)
    return _xgboost(adata, test=test, num_round=num_round)


@_xgboost_method(
    method_name="XGBoost (log scran)",
    image="openproblems-r-extras",
)
def xgboost_scran(adata, test: bool = False, num_round: Optional[int] = None):
    adata = log_scran_pooling(adata)
    return _xgboost(adata, test=test, num_round=num_round)
