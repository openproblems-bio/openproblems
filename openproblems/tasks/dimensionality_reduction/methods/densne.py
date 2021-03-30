from ....tools.decorators import method
from ....tools.utils import check_version

import densne


@method(
    method_name="denSNE",
    paper_name="Assessing single-cell transcriptomic variability through density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/hhcho/densvis",
    code_version="8efe0a2",
)
def densne(adata):
    adata.obsm["X_emb"] = densne.run_densne(adata.X, final_dens=False)
    return adata
