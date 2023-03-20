from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp
from typing import Optional


@method(
    method_name="DestVI",
    method_summary=(
        "destVI is a decomposition method that leverages a conditional generative model"
        " of spatial transcriptomics down to the sub-cell-type variation level, which"
        " is then used to decompose the cell-type proportions determining the spatial"
        " organization of a tissue."
    ),
    paper_name=(
        "DestVI identifies continuums of cell types in spatial transcriptomics data"
    ),
    paper_reference="lopez2022destvi",
    paper_year=2022,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-pytorch",
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
        max_epochs_sc = max_epochs_sc or 500
        max_epochs_sp = max_epochs_sp or 10000

    adata_sc, adata = split_sc_and_sp(adata)

    CondSCVI.setup_anndata(adata_sc, labels_key="label")
    sc_model = CondSCVI(adata_sc, weight_obs=False)
    sc_model.train(
        max_epochs=max_epochs_sc,
        early_stopping=True,
        train_size=0.9,
        validation_size=0.1,
        early_stopping_monitor="elbo_validation",
    )
    DestVI.setup_anndata(adata)

    st_model = DestVI.from_rna_model(adata, sc_model)
    st_model.train(
        max_epochs=max_epochs_sp,
        batch_size=min(int(adata.n_obs / 20 + 3), 128),
        plan_kwargs={"min_kl_weight": 3.0, "max_kl_weight": 3},
    )
    adata.obsm["proportions_pred"] = st_model.get_proportions().to_numpy()
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata
