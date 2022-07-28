from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_r_version

import functools

_r_liana = r_function("liana.R", args="sce, test, ...")

_liana_method = functools.partial(
    method,
    paper_name="Comparison of methods and resources for cell-cell "
               "communication inference from single-cell RNA-Seq data",
    paper_url="https://www.nature.com/articles/s41467-022-30755-0",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_liana_method(
    method_name="LIANA",
)
def liana(adata, score_col="aggregate_rank", asc=True, test=False, **kwargs):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers["logcounts"] = adata.layers["log_cpm"]
    del adata.layers["log_cpm"]

    # Run LIANA
    liana_res = _r_liana(adata, test=test, **kwargs)

    # Format results
    liana_res["score"] = liana_res[score_col]
    liana_res.sort_values("score", ascending=asc, inplace=True)
    adata.uns["ccc"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata


_cellphonedb_method = functools.partial(
    method,
    paper_name="CellPhoneDB: inferring cell–cell communication from "
               "combined expression of multi-subunit ligand–receptor complexes",
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
    adata = liana(
        adata, method="cellphonedb", score_col="lr.mean", asc=False, test=test
    )
    # Filter & Re-order
    adata.uns["ccc"]["score"] = adata.uns["ccc"].apply(
        lambda x: _p_filt(x.pvalue, x["lr.mean"]), axis=1
    )
    adata.uns["ccc"].sort_values("score", ascending=False, inplace=True)

    return adata


_connectome_method = functools.partial(
    method,
    paper_name="Computation and visualization of cell–cell signaling "
               "topologies in single-cell systems data using Connectome",
    paper_url="https://www.nature.com/articles/s41598-022-07959-x",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_connectome_method(
    method_name="Connectome",
)
def connectome(adata, test=False):
    return liana(adata, method="connectome", score_col="weight_sc", asc=False)


_logfc_method = functools.partial(
    method,
    paper_name="Comparison of methods and resources for cell-cell "
    "communication inference from single-cell RNA-Seq data",
    paper_url="https://www.nature.com/articles/s41467-022-30755-0",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_logfc_method(
    method_name="Mean log2FC",
)
def logfc(adata, test=False):
    return liana(adata, method="logfc", score_col="logfc_comb", asc=False)


_natmi_method = functools.partial(
    method,
    paper_name="Predicting cell-to-cell communication networks using NATMI",
    paper_url="https://www.nature.com/articles/s41467-020-18873-z",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_natmi_method(
    method_name="NATMI",
)
def natmi(adata, test=False):
    return liana(adata, method="natmi", score_col="edge_specificity", asc=False)


_sca_method = functools.partial(
    method,
    paper_name="SingleCellSignalR: inference of intercellular networks "
               "from single-cell transcriptomics",
    paper_url="https://academic.oup.com/nar/article/48/10/e55/5810485",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_sca_method(
    method_name="SingleCellSignalR",
)
def sca(adata, test=False):
    return liana(adata, method="sca", score_col="LRscore", asc=False)
