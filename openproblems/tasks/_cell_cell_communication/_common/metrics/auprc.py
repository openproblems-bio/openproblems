from .....tools.decorators import metric
from ..utils import join_truth_and_pred


@metric(
    metric_name="Precision-recall AUC",
    metric_summary="Area under the precision-recall curve.",
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
