from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import split_sc_and_sp


@method(
    method_name="Stereoscope",
    paper_name="Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography",  # noqa: E501
    paper_url="https://www.nature.com/articles/s42003-020-01247-y",
    paper_year=2020,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi-tools"),
    image="openproblems-python-extras",
)
def stereoscope_raw(adata, test=False):
    from scvi.external import RNAStereoscope
    from scvi.external import SpatialStereoscope

    adata_sc, adata = split_sc_and_sp(adata)
    RNAStereoscope.setup_anndata(adata_sc, labels_key="label", layer=None)
    sc_model = RNAStereoscope(adata_sc)
    sc_model.train(max_epochs=100)
    SpatialStereoscope.setup_anndata(adata, layer=None)

    stereo = SpatialStereoscope.from_rna_model(adata, sc_model)
    stereo.train(max_epochs=1000)
    adata.obsm["proportions_pred"] = stereo.get_proportions()
    return adata
