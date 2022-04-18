from ....tools.decorators import method
from ....tools.utils import check_version


@method(
    method_name="DestVI",
    paper_name="Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of inflammation",  # noqa: E501
    paper_url="https://www.biorxiv.org/content/10.1101/2021.05.10.443517v1",
    paper_year=2022,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi-tools"),
    image="openproblems-python-extras",
)
def destvi_raw(adata, test=False):
    from scvi.model import CondSCVI
    from scvi.model import DestVI

    adata_sc = adata.uns["sc_reference"].copy()
    CondSCVI.setup_anndata(adata_sc, labels_key="label", layer=None)
    sc_model = CondSCVI(adata_sc, weight_obs=False)
    sc_model.train()
    DestVI.setup_anndata(adata, layer=None)

    st_model = DestVI.from_rna_model(adata, sc_model)
    st_model.train(max_epochs=2500)
    adata.obsm["proportions_pred"] = st_model.get_proportions()
    return adata
