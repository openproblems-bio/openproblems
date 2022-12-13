from .....tools.decorators import metric
from ..utils import join_truth_and_pred

import numpy as np


def _sigmoid_transform(x):
    return 1 - 1 / (1 + x / 2)


@metric(metric_name="Odds Ratio", maximize=True)
def odds_ratio(adata, top_prop=0.05):
    # Join benchmark (assumed truth) and ccc results
    # Get /w ccc_target and a response [0, 1] column
    gt = join_truth_and_pred(adata)
    gt = gt.sort_values("score", ascending=False)
    top_n = int(adata.uns["ccc_target"].shape[0] * top_prop)

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

    numerator = tp * tn
    denominator = fp * fn
    if denominator == 0:
        if numerator == 0:
            # undefined
            return np.nan
        else:
            # perfect score
            oddsratio = np.inf
    else:
        oddsratio = numerator / denominator

    return _sigmoid_transform(oddsratio)
