from .....tools.decorators import metric

"""
For the bio-conservation score, the ASW was computed on cell identity labels and
scaled to a value between 0 and 1 using the equation:
celltypeASW=(ASWC+1)/2,

where C denotes the set of all cell identity labels."""


@metric(
    metric_name="Silhouette",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def silhouette(adata):
    from scib.metrics import silhouette

    return silhouette(adata, group_key="labels", embed="X_emb")
