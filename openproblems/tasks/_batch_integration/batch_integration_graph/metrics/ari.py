from .....tools.decorators import metric

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
    image="openproblems-python-batch-integration",
)
def ari(adata):
    from scib.metrics import ari
    from scib.metrics.clustering import opt_louvain

    opt_louvain(
        adata,
        label_key="labels",
        cluster_key="cluster",
        plot=False,
        inplace=True,
        force=True,
    )
    return ari(adata, group1="cluster", group2="labels")
