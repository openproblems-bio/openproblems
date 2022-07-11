from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp


@method(
    method_name="Stereoscope",
    paper_name="Single-cell and spatial transcriptomics enables "
    "probabilistic inference of cell type topography",
    paper_url="https://www.nature.com/articles/s42003-020-01247-y",
    paper_year=2020,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-extras",
)
def stereoscope_raw(adata, test=False, max_epochs=None):
    from scvi.external import RNAStereoscope
    from scvi.external import SpatialStereoscope

    if test:
        max_epochs_rna = max_epochs or 10
        max_epochs_spatial = max_epochs or 10
    else:
        max_epochs_rna = max_epochs or 100
        max_epochs_spatial = max_epochs or 1000

    adata_sc, adata = split_sc_and_sp(adata)
    RNAStereoscope.setup_anndata(adata_sc, labels_key="label", layer=None)
    sc_model = RNAStereoscope(adata_sc)
    sc_model.train(max_epochs=max_epochs_rna)
    SpatialStereoscope.setup_anndata(adata, layer=None)

    stereo = SpatialStereoscope.from_rna_model(adata, sc_model)
    stereo.train(max_epochs=max_epochs_spatial)
    adata.obsm["proportions_pred"] = stereo.get_proportions()

    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata
