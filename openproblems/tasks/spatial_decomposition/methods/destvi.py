from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp
from typing import Optional


@method(
    method_name="DestVI",
    paper_name="DestVI identifies continuums of cell types in spatial "
    "transcriptomics data",
    paper_reference="lopez2022destvi",
    paper_year=2022,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-extras",
)
def destvi(
    adata,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_sp: Optional[int] = None,
):
    from scvi.model import CondSCVI
    from scvi.model import DestVI

    if test:
        max_epochs_sp = max_epochs_sp or 10
        max_epochs_sc = max_epochs_sc or 10
    else:  # pragma: nocover
        max_epochs_sc = max_epochs_sc or 300
        max_epochs_sp = max_epochs_sp or 2500

    adata_sc, adata = split_sc_and_sp(adata)

    CondSCVI.setup_anndata(adata_sc, labels_key="label")
    sc_model = CondSCVI(adata_sc, weight_obs=False)
    sc_model.train(
        max_epochs=max_epochs_sc,
        early_stopping=True,
        early_stopping_monitor="reconstruction_loss_train",
    )
    DestVI.setup_anndata(adata)

    st_model = DestVI.from_rna_model(adata, sc_model)
    st_model.train(
        max_epochs=max_epochs_sp,
        early_stopping=True,
        early_stopping_monitor="reconstruction_loss_train",
    )
    adata.obsm["proportions_pred"] = st_model.get_proportions().to_numpy()
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata
