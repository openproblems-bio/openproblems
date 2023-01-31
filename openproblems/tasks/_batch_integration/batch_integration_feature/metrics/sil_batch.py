from .....tools.decorators import metric

"""
We consider the absolute silhouette width, s(i), on
batch labels per cell i. Here, 0 indicates that batches are well mixed, and any
deviation from 0 indicates a batch effect:
𝑠batch(𝑖)=|𝑠(𝑖)|.

To ensure higher scores indicate better batch mixing, these scores are scaled by
subtracting them from 1. As we expect batches to integrate within cell identity
clusters, we compute the batchASWj score for each cell label j separately,
using the equation:
batchASW𝑗=1|𝐶𝑗|∑𝑖∈𝐶𝑗1−𝑠batch(𝑖),

where Cj is the set of cells with the cell label j and |Cj| denotes the number of cells
in that set.

To obtain the final batchASW score, the label-specific batchASWj scores are averaged:
batchASW=1|𝑀|∑𝑗∈𝑀batchASW𝑗.

Here, M is the set of unique cell labels."""


@metric(
    metric_name="Batch ASW",
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-pytorch",
)
def silhouette_batch(adata):
    from ...batch_integration_embed.metrics.sil_batch import (
        silhouette_batch as embed_metric
    )
    from scanpy.tl import pca
    adata.obsm["X_emb"] = pca(adata.X)
    return embed_metric(adata)
