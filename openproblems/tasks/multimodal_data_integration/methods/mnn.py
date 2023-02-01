from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_r_version

import functools

_mnn = r_function("mnn.R")

_mnn_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name=(
        "Batch effects in single-cell RNA-sequencing data are corrected by "
        "matching mutual nearest neighbors"
    ),
    paper_reference="haghverdi2018batch",
    paper_year=2018,
    code_url="https://github.com/LTLA/batchelor",
    image="openproblems-r-extras",
)


@_mnn_method(
    method_name="Mutual Nearest Neighbors (log CPM)",
)
def mnn_log_cpm(adata, test=False):
    adata = log_cpm(adata)
    adata = log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    adata = _mnn(adata)
    adata.uns["method_code_version"] = check_r_version("batchelor")
    return adata


@_mnn_method(
    method_name="Mutual Nearest Neighbors (log scran)",
)
def mnn_log_scran_pooling(adata, test=False):
    adata = log_scran_pooling(adata)
    adata = log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    adata = _mnn(adata)
    adata.uns["method_code_version"] = check_r_version("batchelor")
    return adata
