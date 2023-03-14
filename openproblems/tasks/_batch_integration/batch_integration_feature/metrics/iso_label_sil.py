from .....tools.decorators import metric
from ...batch_integration_embed import metrics as embed_metrics
from .utils import feature_to_embedding

"""
Isolated cell labels are defined as the labels present in the least number
of batches in the integration task. The score evaluates how well these isolated labels
separate from other cell identities.

The isolated label ASW score is obtained by computing the
ASW of isolated versus non-isolated labels on the PCA embedding (ASW metric above) and
scaling this score to be between 0 and 1. The final score for each metric version
consists of the mean isolated score of all isolated labels.
"""


@metric(**embed_metrics.isolated_labels_sil.metadata)
def isolated_labels_sil(adata):
    return embed_metrics.isolated_labels_sil(feature_to_embedding(adata))
