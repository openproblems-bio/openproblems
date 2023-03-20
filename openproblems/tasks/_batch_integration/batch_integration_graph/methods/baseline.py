from .....tools.decorators import baseline_method
from .....tools.utils import check_version
from ..._common.methods.baseline import _random_embedding

import scanpy as sc


@baseline_method(
    method_name="Random Graph by Celltype",
    method_summary=(
        "Cells are embedded as a one-hot encoding of celltype labels. A graph is then"
        " built on this embedding"
    ),
)
def celltype_random_graph(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"])
    sc.pp.neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
