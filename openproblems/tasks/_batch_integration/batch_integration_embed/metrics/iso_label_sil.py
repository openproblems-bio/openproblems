from .....tools.decorators import metric

"""
Isolated cell labels are defined as the labels present in the least number
of batches in the integration task. The score evaluates how well these isolated labels
separate from other cell identities.

The isolated label ASW score is obtained by computing the
ASW of isolated versus non-isolated labels on the PCA embedding (ASW metric above) and
scaling this score to be between 0 and 1. The final score for each metric version
consists of the mean isolated score of all isolated labels.
"""


@metric(
    metric_name="Isolated label Silhouette",
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-pytorch",
)
def isolated_labels_sil(adata):
    from scib.metrics import isolated_labels

    return isolated_labels(
        adata,
        label_key="labels",
        batch_key="batch",
        embed="X_emb",
        cluster=False,
        verbose=False,
    )
