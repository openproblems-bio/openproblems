
from ....tools.decorators import metric
from ..checks import check_method


@metric(metric_name="Pearson's chi squared test", maximize=True)
def chi_squared(adata, method):

    """
    Perform Pearson's chi squared test to estimate the goodness of fit of dimensional reduction.
    """
    
    if check_method(adata) == True:

        from scipy.stats import chisquare

        adata.obsp["chi_squared"] = chisquare(adata.obsm["dimensional_reduction_method"], f_exp=adata.X.astype(array))

        ### need to add a summation method to this

        return float(adata.uns["chi_squared_sum"])