from ....tools.decorators import method


@method(
    method_name="densMAP",
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/lmcinnes/umap",
    code_version="8efe0a2",
    image="openproblems-python-extras",
)
def densmap(adata):
    import umap

    adata.obsm["X_emb"] = umap.UMAP(densmap=True, random_state=42).fit_transform(
        adata.X
    )
    return adata  # return with embedding
