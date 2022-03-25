from .....tools.decorators import metric

"""NMI compares the overlap of two clusterings.
We used NMI to compare the cell-type labels with Louvain clusters computed on the integrated dataset.
The overlap was scaled using the mean of the entropy terms for cell-type and cluster labels.
Thus, NMI scores of 0 or 1 correspond to uncorrelated clustering or a perfect match, respectively.
We performed optimized Louvain clustering for this metric to obtain the best match between clusters and labels.
Louvain clustering was performed at a resolution range of 0.1 to 2 in steps of 0.1, and the clustering output
with the highest NMI with the label set was used. We used the scikit-learn27 (v.0.22.1) implementation of NMI.
"""


@metric(
    metric_name="NMI",
    maximize=True,
    image="openproblems-python-batch-integration",
)
def nmi(adata):
    from scib.metrics.clustering import opt_louvain  # isort:skip
    from scib.metrics import nmi

    opt_louvain(
        adata,
        label_key="labels",
        cluster_key="cluster",
        plot=False,
        inplace=True,
        force=True,
    )
    return nmi(adata, group1="cluster", group2="labels")
