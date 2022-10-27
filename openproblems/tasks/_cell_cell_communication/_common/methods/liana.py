from .....tools.conversion import r_function
from .....tools.decorators import method
from .....tools.normalize import log_cpm
from .....tools.utils import check_r_version
from ..utils import aggregate_method_scores
from ..utils import ligand_receptor_resource

import functools


# Helper function to filter according to permutation p-values
def _p_filt(x, y):
    if x <= 0.05:
        return y
    else:
        return 0


_r_liana = r_function(
    "liana.R", args="sce, op_resource, min_expression_prop, idents_col, test, ..."
)

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
    method_name="LIANA's RobustRankAggregate Rank",
)
def liana(
    adata,
    score_col="aggregate_rank",
    min_expression_prop=0.1,
    test=False,
    **kwargs,
):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers["logcounts"] = adata.layers["log_cpm"]
    del adata.layers["log_cpm"]

    # Run LIANA
    liana_res = _r_liana(
        adata,
        op_resource=ligand_receptor_resource(adata.uns["target_organism"]),
        min_expression_prop=min_expression_prop,
        idents_col="label",
        test=test,
        **kwargs,
    )

    # Format results
    liana_res["score"] = liana_res[score_col]
    adata.uns["ccc_pred"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata


@_liana_method(
    method_name="RobustRankAggregate MAX",
)
def liana_max(adata, test=False):
    adata = liana(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_liana_method(
    method_name="RobustRankAggregate SUM",
)
def liana_sum(adata, test=False):
    adata = liana(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

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


@_cellphonedb_method(
    method_name="CellPhoneDB",
)
def cellphonedb(adata, test=False):
    adata = liana(
        adata,
        method="cellphonedb",
        score_col="lr.mean",
        test=test,
        complex_policy="min",
    )
    # Filter & Re-order
    adata.uns["ccc_pred"]["score"] = adata.uns["ccc_pred"].apply(
        lambda x: _p_filt(x.pvalue, x["lr.mean"]), axis=1
    )

    return adata


@_cellphonedb_method(
    method_name="CellPhoneDB MAX",
)
def cellphonedb_max(adata, test=False):
    adata = cellphonedb(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_cellphonedb_method(
    method_name="CellPhoneDB SUM",
)
def cellphonedb_sum(adata, test=False):
    adata = cellphonedb(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

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
    return liana(adata, method="connectome", score_col="weight_sc", test=test)


@_connectome_method(
    method_name="Connectome MAX",
)
def connectome_max(adata, test=False):
    adata = connectome(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_connectome_method(
    method_name="Connectome SUM",
)
def connectome_sum(adata, test=False):
    adata = connectome(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


@_liana_method(
    method_name="Mean log2FC",
)
def logfc(adata, test=False):
    return liana(adata, method="logfc", score_col="logfc_comb", test=test)


@_connectome_method(
    method_name="Log2FC MAX",
)
def logfc_max(adata, test=False):
    adata = logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_connectome_method(
    method_name="Log2FC SUM",
)
def logfc_sum(adata, test=False):
    adata = logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


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
    return liana(adata, method="natmi", score_col="edge_specificity", test=test)


@_natmi_method(
    method_name="NATMI MAX",
)
def natmi_max(adata, test=False):
    adata = natmi(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_natmi_method(
    method_name="NATMI SUM",
)
def natmi_sum(adata, test=False):
    adata = natmi(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


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
    return liana(adata, method="sca", score_col="LRscore", test=test)


@_sca_method(
    method_name="SingleCellSignalR MAX",
)
def sca_max(adata, test=False):
    adata = sca(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_sca_method(
    method_name="SingleCellSignalR SUM",
)
def sca_sum(adata, test=False):
    adata = sca(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata
