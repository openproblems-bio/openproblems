from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import split_sc_and_sp


@method(
    method_name="DestVI",
    paper_name="Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of inflammation",  # noqa: E501
    paper_url="https://www.biorxiv.org/content/10.1101/2021.05.10.443517v1",
    paper_year=2022,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-extras",
)
def destvi_raw(adata, test=False, max_epochs=None):
    from scvi.model import CondSCVI
    from scvi.model import DestVI

    if test:
        max_epochs = max_epochs or 10
    else:
        max_epochs = max_epochs or 2500

    adata_sc, adata_sp = split_sc_and_sp(adata)
    CondSCVI.setup_anndata(adata_sc, labels_key="label", layer=None)
    sc_model = CondSCVI(adata_sc, weight_obs=False)
    sc_model.train()
    DestVI.setup_anndata(adata_sp, layer=None)

    st_model = DestVI.from_rna_model(adata_sp, sc_model)
    st_model.train(max_epochs=max_epochs)
    adata_sp.obsm["proportions_pred"] = st_model.get_proportions()
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata_sp
