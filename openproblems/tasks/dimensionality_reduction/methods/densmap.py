from ....tools.decorators import method

import densmap as dens


@method(
    method_name="densMAP",
    paper_name="Assessing single-cell transcriptomic variability through"
    " density-preserving data visualization",
    paper_url="https://www.nature.com/articles/s41587-020-00801-7",
    paper_year=2021,
    code_url="https://github.com/hhcho/densvis",
    code_version="8efe0a2",
    image="openproblems-python-extras"
)
def densmap(adata):
    adata.obsm["X_emb"] = dens.densMAP(final_dens=False).fit_transform(adata.X)
    return adata
