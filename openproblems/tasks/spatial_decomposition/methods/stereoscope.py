from ....tools.decorators import method
from ....tools.utils import check_version


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
    from scvi.data import setup_anndata
    from scvi.external.stereoscope import RNAStereoscope
    from scvi.external.stereoscope import SpatialStereoscope

    adata_sc = adata.uns["sc_reference"].copy()
    setup_anndata(adata_sc, labels_key="label", layer=None)
    sc_model = RNAStereoscope(adata_sc)
    sc_model.train()
    setup_anndata(adata, layer=None)

    stereo = SpatialStereoscope.from_rna_model(adata, sc_model)
    stereo.train()
    adata.obsm["proportions_pred"] = stereo.get_proportions()
    return adata
