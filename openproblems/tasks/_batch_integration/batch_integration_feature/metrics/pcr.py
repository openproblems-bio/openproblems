from .....tools.decorators import metric

"""
Principal component regression, derived from PCA, has previously been used to quantify
batch removal. Briefly, the R2 was calculated from a linear regression of the
covariate of interest (for example, the batch variable B) onto each principal component.
The variance contribution of the batch effect per principal component was then
calculated as the product of the variance explained by the ith principal component (PC)
and the corresponding R2(PCi|B). The sum across all variance contributions by the batch
effects in all principal components gives the total variance explained by the batch
variable as follows:
Var(𝐶|𝐵)=∑𝑖=1𝐺Var(𝐶|PC𝑖)×𝑅2(PC𝑖|𝐵),

where Var(C|PCi) is the variance of the data matrix C explained by the ith principal
component."""


@metric(
    metric_name="PC Regression",
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-pytorch",
)
def pcr(adata):
    from ...batch_integration_embed.metrics.pcr import (
        pcr as embed_metric
    )
    from scanpy.tl import pca
    adata.obsm["X_emb"] = pca(adata.X)
    return embed_metric(adata)
