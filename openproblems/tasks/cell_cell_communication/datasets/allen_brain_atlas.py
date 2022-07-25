from ....data.allen_brain_atlas import load_mouse_brain_atlas
from ....tools.decorators import dataset


@dataset(
    "~15k Murine Brain cells" 
    "Tasic, B. et al., 2016. Nat. Neurosci. 19, 335–346 .",
    data_url=load_mouse_brain_atlas.metadata["data_url"],
    dataset_summary="TODO",
    image="openproblems-r-extras",
)
def mouse_brain_atlas(test=False):
    return load_mouse_brain_atlas(test=test)
