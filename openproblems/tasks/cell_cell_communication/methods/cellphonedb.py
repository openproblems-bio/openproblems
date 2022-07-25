from ....tools.decorators import method
from .liana import liana

import functools

_cellphonedb_method = functools.partial(
    method,
    paper_name="CellPhoneDB: inferring cell–cell communication from "
    "combined expression of multi-subunit ligand–receptor complexes.",
    paper_url="https://www.nature.com/articles/s41596-020-0292-x",
    paper_year=2020,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


# Helper function to filter according to permutation p-values
def _p_filt(x, y):
    if x <= 0.05:
        return y
    else:
        return 0


@_cellphonedb_method(
    method_name="CellPhoneDB",
)
def cellphonedb(adata, test=False):

    adata = liana(adata,
                  method="cellphonedb",
                  score_col="lr.mean",
                  asc=False,
                  test=test)
    # Filter & Re-order
    adata.uns["ccc"]["score"] = adata.uns["ccc"].apply(
        lambda x: _p_filt(x.pvalue, x["lr.mean"]), axis=1
    )
    adata.uns["ccc"].sort_values("score", ascending=False, inplace=True)

    return adata
