from .....tools.decorators import metric
from ..utils import join_truth_and_pred
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve


@metric(metric_name="Precision-recall AUC", maximize=True)
def auprc(adata):
    gt = join_truth_and_pred(adata)
    precision, recall, _ = precision_recall_curve(
        gt["response"], gt["score"], pos_label=1
    )
    return auc(recall, precision)
