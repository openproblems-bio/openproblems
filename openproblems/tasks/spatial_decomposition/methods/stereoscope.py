from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp


@method(
    method_name="Stereoscope",
    paper_name="Single-cell and spatial transcriptomics enables probabilistic "
    "inference of cell type topography",
    paper_reference="andersson2020single",
    paper_year=2020,
    code_url="https://github.com/scverse/scvi-tools",
    image="openproblems-python-pytorch",
)
def stereoscope(adata, test=False, max_epochs_sc=None, max_epochs_sp=None):
    from scvi.external import RNAStereoscope
    from scvi.external import SpatialStereoscope

    if test:
        max_epochs_sp = max_epochs_sp or 10
        max_epochs_sc = max_epochs_sc or 10
    else:  # pragma: nocover
        max_epochs_sc = max_epochs_sc or 100
        max_epochs_sp = max_epochs_sp or 1000

    adata_sc, adata = split_sc_and_sp(adata)

    RNAStereoscope.setup_anndata(adata_sc, labels_key="label")
    sc_model = RNAStereoscope(adata_sc)
    sc_model.train(
        max_epochs=max_epochs_sc,
        early_stopping=True,
        early_stopping_monitor="elbo_train",
    )

    SpatialStereoscope.setup_anndata(adata)
    st_model = SpatialStereoscope.from_rna_model(adata, sc_model)
    st_model.train(
        max_epochs=max_epochs_sp,
        early_stopping=True,
        early_stopping_monitor="elbo_train",
    )
    adata.obsm["proportions_pred"] = st_model.get_proportions().to_numpy()
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata
