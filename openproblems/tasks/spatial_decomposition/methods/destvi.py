from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp


@method(
    method_name="DestVI",
    paper_name="DestVI identifies continuums of cell types in spatial transcriptomics data",  # noqa: E501
    paper_url="https://doi.org/10.1038/s41587-022-01272-8",
    paper_year=2022,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-extras",
)
def destvi(adata, test=False, max_epochs_sc=None, max_epochs_sp=None):
    from scvi.model import CondSCVI
    from scvi.model import DestVI

    if test:
        max_epochs_sp = max_epochs_sp or 10
        max_epochs_sc = max_epochs_sc or 10
    else:
        max_epochs_sc = max_epochs_sc or 300
        max_epochs_sp = max_epochs_sp or 2500

    adata_sc, adata_sp = split_sc_and_sp(adata)

    CondSCVI.setup_anndata(adata_sc, labels_key="label")
    sc_model = CondSCVI(adata_sc, weight_obs=False)
    sc_model.train(max_epochs=max_epochs_sc)
    DestVI.setup_anndata(adata_sp)

    st_model = DestVI.from_rna_model(adata_sp, sc_model)
    st_model.train(max_epochs=max_epochs_sp)
    adata_sp.obsm["proportions_pred"] = st_model.get_proportions()
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata_sp
