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
    paper_reference="dimitrov2022comparison",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


def _liana(
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
    method_name="LIANA Rank Aggregate (max)",
)
def liana_max(adata, test=False):
    adata = _liana(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_liana_method(
    method_name="LIANA Rank Aggregate (sum)",
)
def liana_sum(adata, test=False):
    adata = _liana(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_cellphonedb_method = functools.partial(
    method,
    paper_name="CellPhoneDB: inferring cell–cell communication from "
    "combined expression of multi-subunit ligand–receptor complexes",
    paper_reference="efremova2020cellphonedb",
    paper_year=2020,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


def _cellphonedb(adata, test=False):
    adata = _liana(
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
    method_name="CellPhoneDB (max)",
)
def cellphonedb_max(adata, test=False):
    adata = _cellphonedb(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_cellphonedb_method(
    method_name="CellPhoneDB (sum)",
)
def cellphonedb_sum(adata, test=False):
    adata = _cellphonedb(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_connectome_method = functools.partial(
    method,
    paper_name="Computation and visualization of cell–cell signaling "
    "topologies in single-cell systems data using Connectome",
    paper_reference="raredon2022computation",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


def _connectome(adata, test=False):
    return _liana(adata, method="connectome", score_col="weight_sc", test=test)


@_connectome_method(
    method_name="Connectome (max)",
)
def connectome_max(adata, test=False):
    adata = _connectome(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_connectome_method(
    method_name="Connectome (sum)",
)
def connectome_sum(adata, test=False):
    adata = _connectome(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


def _logfc(adata, test=False):
    return _liana(adata, method="logfc", score_col="logfc_comb", test=test)


@_connectome_method(
    method_name="Log2FC (max)",
)
def logfc_max(adata, test=False):
    adata = _logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_connectome_method(
    method_name="Log2FC (sum)",
)
def logfc_sum(adata, test=False):
    adata = _logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_natmi_method = functools.partial(
    method,
    paper_name="Predicting cell-to-cell communication networks using NATMI",
    paper_reference="hou2020predicting",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


def _natmi(adata, test=False):
    return _liana(adata, method="natmi", score_col="edge_specificity", test=test)


@_natmi_method(
    method_name="NATMI (max)",
)
def natmi_max(adata, test=False):
    adata = _natmi(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_natmi_method(
    method_name="NATMI (sum)",
)
def natmi_sum(adata, test=False):
    adata = _natmi(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_sca_method = functools.partial(
    method,
    paper_name="SingleCellSignalR: inference of intercellular networks "
    "from single-cell transcriptomics",
    paper_reference="cabello2020singlecellsignalr",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


def _sca(adata, test=False):
    return _liana(adata, method="sca", score_col="LRscore", test=test)


@_sca_method(
    method_name="SingleCellSignalR (max)",
)
def sca_max(adata, test=False):
    adata = _sca(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_sca_method(
    method_name="SingleCellSignalR (sum)",
)
def sca_sum(adata, test=False):
    adata = _sca(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata
