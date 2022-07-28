from ....data.allen_brain_atlas import load_mouse_brain_atlas
from ....tools.decorators import dataset


@dataset(
    "~15k Murine Brain cells",
    data_url=load_mouse_brain_atlas.metadata["data_url"],
    data_reference=load_mouse_brain_atlas.metadata["data_reference"],
    dataset_summary="A murine brain atlas with inferred spatially-adjacent "
    "cell types as assumed benchmark truth. Adjacent cell types are those "
    "with z-transformed deconvolution proportion correlations >= 1.645, "
    "inferred on murine brain 10x Visium slides (see Dimitrov et al., 2022)",
    image="openproblems-r-extras",
)
def mouse_brain_atlas(test=False):
    return load_mouse_brain_atlas(test=test)
