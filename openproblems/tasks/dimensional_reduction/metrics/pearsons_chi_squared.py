
from ....tools.decorators import metric


@metric(metric_name="Pearson's chi squared test", maximize=True)
def chi_squared(adata, method):

    """
    Perform Pearson's chi squared test to estimate the goodness of fit of dimensional reduction.
    """

    from scipy.stats import chisquare

    adata.obsp["chi_squared"] = chisquare(adata.obsm["dimensional_reduction_method"], f_exp=adata.X.astype(array))

