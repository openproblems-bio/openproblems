from .....tools.decorators import method
from .....tools.utils import check_version
from ..._common.methods.baseline import _random_embedding

import scanpy as sc


@method(
    method_name="Random Graph by Celltype",
    paper_name="Random Graph by Celltype (baseline)",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def celltype_random_graph(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"])
    sc.pp.neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
