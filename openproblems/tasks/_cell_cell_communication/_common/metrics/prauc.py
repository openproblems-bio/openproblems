from .....tools.decorators import metric
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve

import numpy as np


@metric(metric_name="Precision-recall AUC", maximize=True)
def prauc(adata, merge_keys, n_times=100):
    gt = adata.uns["ccc_target"].merge(
        adata.uns["ccc_pred"], on=merge_keys, how="outer"
    )

    gt.loc[gt["response"].isna(), "response"] = 0
    gt.loc[gt["score"].isna(), "score"] = np.nanmin(gt["score"]) - 1

    cpn = np.sum(gt["response"] == 1)  # condition positive num
    cnn = np.sum(gt["response"] == 0)  # condition negative num

    if (cpn != cnn) & (n_times > 1):
        sub_n = np.min((cpn, cnn))
        pr = _auprcs_sub(gt=gt, sub_n=sub_n, n_times=n_times)
    else:
        pr = _get_auc(gt["response"], gt["score"], pos_label=1)
    return pr


def _get_auc(response, score, pos_label):
    precision, recall, _ = precision_recall_curve(response, score, pos_label=pos_label)
    return auc(recall, precision)


def _auprcs_sub(gt, sub_n, n_times, seed=1):
    pos_idx = np.where(gt["response"] == 1)[0]
    neg_idx = np.where(gt["response"] == 0)[0]

    # rng instance
    rng = np.random.default_rng(seed=seed)
    prcs = []

    for x in range(n_times):
        # join /w shuffled balance size indices
        sneg_idx = rng.choice(neg_idx, sub_n, replace=False)
        spos_idex = rng.choice(pos_idx, sub_n, replace=False)
        bal_idx = np.union1d(sneg_idx, spos_idex)

        temp = gt.iloc[bal_idx].copy()
        prcs.append(_get_auc(temp["response"], temp["score"], pos_label=1))
    pr = np.mean(prcs)

    return pr
