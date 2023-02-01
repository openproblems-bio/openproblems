from .....tools.decorators import metric
from ...batch_integration_graph import metrics as graph_metrics

"""
The Rand index compares the overlap of two clusterings;
it considers both correct clustering overlaps while also counting correct
disagreements between two clusterings.
Similar to NMI, we compared the cell-type labels with the NMI-optimized
Louvain clustering computed on the integrated dataset.
The adjustment of the Rand index corrects for randomly correct labels.
An ARI of 0 or 1 corresponds to random labeling or a perfect match,
respectively.
We also used the scikit-learn (v.0.22.1) implementation of the ARI.
"""


@metric(
    metric_name="ARI",
    maximize=True,
    paper_reference="luecken2022benchmarking",
    image="openproblems-r-pytorch",
)
def ari(adata):
    from scanpy.pp import neighbors

    neighbors(adata, use_rep="X_emb")
    return graph_metrics.ari(adata)
