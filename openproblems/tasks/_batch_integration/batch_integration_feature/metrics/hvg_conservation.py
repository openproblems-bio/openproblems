from .....tools.decorators import metric

"""
The HVG conservation score is a proxy for the preservation of
the biological signal. If the data integration method returned
a corrected data matrix, we computed the number of HVGs before
and after correction for each batch via Scanpy’s
highly_variable_genes function (using the ‘cell ranger’ flavor).
If available, we computed 500 HVGs per batch. If fewer than 500
genes were present in the integrated object for a batch,
the number of HVGs was set to half the total genes in that batch.
The overlap coefficient is as follows:
overlap(𝑋,𝑌)=|𝑋∩𝑌|/min(|𝑋|,|𝑌|),

where X and Y denote the fraction of preserved informative genes.
The overall HVG score is the mean of the per-batch HVG overlap
coefficients.
"""


@metric(
    metric_name="HVG conservation",
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-python-batch-integration",
)
def hvg_conservation(adata):
    from scib.metrics import hvg_overlap

    adata_unint = adata.copy()
    adata_unint.X = adata_unint.layers["log_normalized"]

    return hvg_overlap(adata_unint, adata, "batch")
