from ....tools.decorators import method
from ....tools.utils import check_version
from .preprocessing import preprocess_logCPM_1kHVG


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
def densmap_logCPM_1kHVG(adata):
    from umap import UMAP

    preprocess_logCPM_1kHVG(adata)
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
def densmap_pca_logCPM_1kHVG(adata):
    from umap import UMAP

    preprocess_logCPM_1kHVG(adata)
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(
        adata.obsm["X_input"]
    )
    return adata
