from .....tools.decorators import metric

import numpy as np


@metric(metric_name="Odds Ratio", maximize=True)
def odds_ratio(adata, merge_keys, top_n=100):
    # Join benchmark (assumed truth) and ccc results
    # Get /w ccc_target and a response [0, 1] column
    gt = adata.uns["ccc_target"].merge(
        adata.uns["ccc_pred"], on=merge_keys, how="inner"
    )

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

    if tp == top_n:  # best case scenario (all are true predictions)
        tp = top_n - 1
        fp = 1
    if tp == 0:  # worst case scenario (no true predictions)
        tp = 1
        fp = top_n - 1

    oddsratio = (tp / fp) / (fn / tn)

    return oddsratio
