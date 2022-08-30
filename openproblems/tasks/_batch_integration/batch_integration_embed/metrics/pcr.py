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
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def pcr(adata):
    from ._utils import _get_split
    from scib.metrics import pcr_comparison

    return pcr_comparison(*_get_split(adata), "batch", embed="X_emb")
