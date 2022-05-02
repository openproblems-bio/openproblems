from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import scanpy as sc


@method(
    method_name="densMAP (logCPM, 1kHVG)",
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    code_version=check_version("umap-learn"),
    image="openproblems-python-extras",
)
def densmap_logCPM_1kHVG(adata, test: bool = False):
    from umap import UMAP

    adata = log_cpm_hvg(adata)
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(adata.X)
    return adata


@method(
    method_name="densMAP PCA (logCPM, 1kHVG)",
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    code_version=check_version("umap-learn"),
    image="openproblems-python-extras",
)
def densmap_pca_logCPM_1kHVG(adata, test: bool = False):
    from umap import UMAP

    adata = log_cpm_hvg(adata)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(
        adata.obsm["X_pca"]
    )
    return adata
