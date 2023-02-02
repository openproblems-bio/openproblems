from .....tools.decorators import metric
from ..utils import join_truth_and_pred


@metric(
    metric_name="Precision-recall AUC",
    metric_summary=(
        "a single number _[0-1]_ that summarizes the area under the curve where x is"
        " the recall and y is the precision."
    ),
    paper_reference="davis2006prauc",
    maximize=True,
)
def auprc(adata):
    from sklearn.metrics import auc
    from sklearn.metrics import precision_recall_curve

    gt = join_truth_and_pred(adata)
    precision, recall, _ = precision_recall_curve(
        gt["response"], gt["score"], pos_label=1
    )
    return auc(recall, precision)
