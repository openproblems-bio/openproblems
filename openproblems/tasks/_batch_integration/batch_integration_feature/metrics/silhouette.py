from .....tools.decorators import metric
from ...batch_integration_embed import metrics as embed_metrics
from .utils import feature_to_embedding

"""
For the bio-conservation score, the ASW was computed on cell identity labels and
scaled to a value between 0 and 1 using the equation:
celltypeASW=(ASW_C+1)/2,

where C denotes the set of all cell identity labels.
For information about the batch silhouette score, check sil_batch."""


@metric(**embed_metrics.silhouette.metadata)
def silhouette(adata):
    return embed_metrics.silhouette(feature_to_embedding(adata))
