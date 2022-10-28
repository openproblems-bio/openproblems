from .....data.allen_brain_atlas import load_mouse_brain_atlas
from .....tools.decorators import dataset
from ..._common.utils import ligand_receptor_resource


@dataset(
    "Mouse brain atlas",
    data_url=load_mouse_brain_atlas.metadata["data_url"],
    data_reference=load_mouse_brain_atlas.metadata["data_reference"],
    dataset_summary="A murine brain atlas with adjacent cell types as assumed "
                    "benchmark truth, inferred from deconvolution proportion "
                    "correlations using matching 10x Visium slides "
                    "(see Dimitrov et al., 2022)."
                    " 14249 cells x 34617 features with 23 cell type labels.",
    image="openproblems-r-extras",
)
def mouse_brain_atlas(test=False):
    adata = load_mouse_brain_atlas(test=test)
    adata.uns["merge_keys"] = ["source", "target"]
    adata.uns["ligand_receptor_resource"] = ligand_receptor_resource(
        adata.uns["target_organism"]
    )
    return adata
