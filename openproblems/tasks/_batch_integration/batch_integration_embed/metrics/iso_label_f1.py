from .....tools.decorators import metric
from ...batch_integration_graph import metrics as graph_metrics

"""
We developed two isolated label scores to evaluate how well the data integration methods
dealt with cell identity labels shared by few batches. Specifically, we identified
isolated cell labels as the labels present in the least number of batches in the
integration task.
The score evaluates how well these isolated labels separate from other cell identities.
We implemented the isolated label metric in two versions:
(1) the best clustering of the isolated label (F1 score) and
(2) the global ASW of the isolated label. For the cluster-based score,
we first optimize the cluster assignment of the isolated label using the F1 score
across louvain clustering resolutions ranging from 0.1 to 2 in resolution steps of 0.1.
The optimal F1 score for the isolated label is then used as the metric score.
The F1 score is a weighted mean of precision and recall given by the equation:
𝐹1=2×(precision×recall)/(precision+recall).

It returns a value between 0 and 1,
where 1 shows that all of the isolated label cells and no others are captured in
the cluster. For the isolated label ASW score, we compute the ASW of isolated
versus nonisolated labels on the PCA embedding (ASW metric above) and scale this
score to be between 0 and 1. The final score for each metric version consists of
the mean isolated score of all isolated labels.
"""


@metric(**graph_metrics.isolated_labels_f1.metadata)
def isolated_labels_f1(adata):
    from scanpy.pp import neighbors

    neighbors(adata, use_rep="X_emb")
    return graph_metrics.isolated_labels_f1(adata)
