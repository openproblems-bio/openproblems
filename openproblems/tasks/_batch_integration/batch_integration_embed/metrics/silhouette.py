from .....tools.decorators import metric

"""
For the bio-conservation score, the ASW was computed on cell identity labels and
scaled to a value between 0 and 1 using the equation:
celltypeASW=(ASW_C+1)/2,

where C denotes the set of all cell identity labels.
For information about the batch silhouette score, check sil_batch."""


@metric(
    metric_name="Silhouette",
    metric_summary=(
        "The absolute silhouette with is computed on cell identity labels, measuring"
        " their compactness."
    ),
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-pytorch",
)
def silhouette(adata):
    from scib.metrics import silhouette

    return silhouette(adata, group_key="labels", embed="X_emb")
