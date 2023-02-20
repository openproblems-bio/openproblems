from .....tools.conversion import r_function
from .....tools.decorators import method
from .....tools.normalize import log_cp10k
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
    "liana.R",
    args="sce, op_resource, min_expression_prop, idents_col, test, aggregate_how, ...",
)

_liana_method = functools.partial(
    method,
    method_summary=(
        "RobustRankAggregate generates a consensus rank of all methods implemented in"
        " LIANA providing either specificity or magnitude scores."
    ),
    paper_name=(
        "Comparison of methods and resources for cell-cell communication inference from"
        " single-cell RNA-Seq data"
    ),
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
    aggregate_how=None,
    **kwargs,
):
    # log-normalize
    adata = log_cp10k(adata)
    adata.layers["logcounts"] = adata.layers["log_cp10k"]
    del adata.layers["log_cp10k"]

    # Run LIANA
    liana_res = _r_liana(
        adata,
        op_resource=ligand_receptor_resource(adata.uns["target_organism"]),
        min_expression_prop=min_expression_prop,
        idents_col="label",
        test=test,
        aggregate_how=aggregate_how,
        **kwargs,
    )

    # Format results
    liana_res["score"] = liana_res[score_col]
    adata.uns["ccc_pred"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata


@_liana_method(
    method_name="Specificity Rank Aggregate (max)",
)
def specificity_max(adata, test=False):
    adata = _liana(adata, test=test, aggregate_how="specificity")
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_liana_method(
    method_name="Specificity Rank Aggregate (sum)",
)
def specificity_sum(adata, test=False):
    adata = _liana(adata, test=test, aggregate_how="specificity")
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


@_liana_method(
    method_name="Magnitude Rank Aggregate (max)",
)
def magnitude_max(adata, test=False):
    adata = _liana(adata, test=test, aggregate_how="magnitude")
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_liana_method(
    method_name="Magnitude Rank Aggregate (sum)",
)
def magnitude_sum(adata, test=False):
    adata = _liana(adata, test=test, aggregate_how="magnitude")
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_cellphonedb_method = functools.partial(
    _liana_method,
    method_summary=(
        "CellPhoneDBv2 calculates a mean of ligand-receptor expression as a measure of"
        " interaction magnitude, along with a permutation-based p-value as a measure of"
        " specificity. Here, we use the former to prioritize interactions, subsequent"
        " to filtering according to p-value less than 0.05."
    ),
    paper_name=(
        "CellPhoneDB: inferring cell–cell communication from combined expression of"
        " multi-subunit ligand–receptor complexes"
    ),
    paper_reference="efremova2020cellphonedb",
    paper_year=2020,
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
    _liana_method,
    method_summary=(
        "Connectome uses the product of ligand-receptor expression as a measure of"
        " magnitude, and the average of the z-transformed expression of ligand and"
        " receptor as a measure of specificity."
    ),
    paper_name=(
        "Computation and visualization of cell–cell signaling topologies in single-cell"
        " systems data using Connectome"
    ),
    paper_reference="raredon2022computation",
    paper_year=2022,
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


_logfc_method = functools.partial(
    _liana_method,
    method_summary=(
        "logFC (implemented in LIANA and inspired by iTALK) combines both expression"
        " and magnitude, and represents the average of one-versus-the-rest log2-fold"
        " change of ligand and receptor expression per cell type."
    ),
)


def _logfc(adata, test=False):
    return _liana(adata, method="logfc", score_col="logfc_comb", test=test)


@_logfc_method(
    method_name="Log2FC (max)",
)
def logfc_max(adata, test=False):
    adata = _logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="max")

    return adata


@_logfc_method(
    method_name="Log2FC (sum)",
)
def logfc_sum(adata, test=False):
    adata = _logfc(adata, test=test)
    adata.uns["ccc_pred"] = aggregate_method_scores(adata, how="sum")

    return adata


_natmi_method = functools.partial(
    _liana_method,
    method_summary=(
        "NATMI uses the product of ligand-receptor expression as a measure of"
        " magnitude. As a measure of specificity, NATMI proposes $specificity.edge ="
        r" \frac{l}{l_s} \cdot \frac{r}{r_s}$; where $l$ and $r$ represent the average"
        " expression of ligand and receptor per cell type, and $l_s$ and $r_s$"
        " represent the sums of the average ligand and receptor expression across all"
        " cell types. We use its specificity measure, as recommended by the authors for"
        " single-context predictions."
    ),
    paper_name="Predicting cell-to-cell communication networks using NATMI",
    paper_reference="hou2020predicting",
    paper_year=2021,
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
    _liana_method,
    method_summary=(
        "SingleCellSignalR provides a magnitude score as $LRscore ="
        r" \frac{\sqrt{lr}}{\mu+\sqrt{lr}}$; where $l$ and $r$ are the average ligand"
        r" and receptor expression per cell type, and $\mu$ is the mean of the"
        " expression matrix."
    ),
    paper_name=(
        "SingleCellSignalR: inference of intercellular networks from single-cell"
        " transcriptomics"
    ),
    paper_reference="cabello2020singlecellsignalr",
    paper_year=2021,
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
