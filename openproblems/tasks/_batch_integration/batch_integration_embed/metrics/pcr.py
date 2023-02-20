from .....tools.decorators import metric
from .utils import get_split

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
    metric_summary=(
        "This compares the explained variance by batch before and after integration. It"
        " returns a score between 0 and 1 (scaled=True) with 0 if the variance"
        " contribution hasn’t changed. The larger the score, the more different the"
        " variance contributions are before and after integration."
    ),
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-pytorch",
)
def pcr(adata):
    from scib.metrics import pcr_comparison

    return pcr_comparison(*get_split(adata), "batch", embed="X_emb")
