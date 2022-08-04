from ....tools.decorators import metric

import numpy as np
import scipy.stats as stats


@metric(metric_name="Odds Ratio", maximize=True, image="openproblems-r-extras")
def odds_ratio(adata, top_n=100):
    # Join benchmark (assumed truth) and ccc results
    # Get /w ccc_target and a response [0, 1] column
    gt = adata.uns["ccc_target"].merge(adata.uns["ccc_pred"], how="right")
    gt = gt[gt["response"].notna()]

    # assign the top rank interactions to 1
    a = np.zeros(len(gt["score"]))
    a[0:top_n] = 1
    gt.loc[:, ["top_n"]] = a

    # Shape to contingency table
    table = np.array(gt.pivot_table(index=["top_n", "response"], aggfunc="size"))

    # if positive or negative class is not in top_n
    if table.shape != (4,):
        return 1

    # Fisher ET
    oddsratio, _ = stats.fisher_exact(table.reshape(2, 2))

    return oddsratio
