from .....tools.decorators import metric
from ..utils import join_truth_and_pred

import numpy as np


@metric(metric_name="Odds Ratio", maximize=True)
def odds_ratio(adata, top_prop=0.05):
    # Join benchmark (assumed truth) and ccc results
    # Get /w ccc_target and a response [0, 1] column
    gt = join_truth_and_pred(adata)
    gt = gt.sort_values("score", ascending=False)
    top_n = np.int(adata.uns["ccc_target"].shape[0] * top_prop)

    # assign the top rank interactions to 1
    a = np.zeros(len(gt["score"]))
    a[0:top_n] = 1
    gt.loc[:, ["top_n"]] = a

    top = gt[gt["top_n"] == 1]
    tp = np.sum(top.response == 1)
    fp = np.sum(top.response == 0)

    bot = gt[gt["top_n"] == 0]
    fn = np.sum(bot.response == 1)
    tn = np.sum(bot.response == 0)

    # GT dependent
    if np.sum(bot.response) == 0:  # comparison is impossible
        return np.inf

    if tp == top_n:  # best case scenario (all are true predictions)
        tp = top_n - 1
        fp = 1
    if tp == 0:  # worst case scenario (no true predictions)
        tp = 1
        fp = top_n - 1

    oddsratio = (tp / fp) / (fn / tn)

    return oddsratio
