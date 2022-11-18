from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset
from . import utils


@dataset(
    dataset_name="Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreatic islet scRNA-seq data from 6 datasets "
    "across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, "
    "and SMARTER-seq). Here we just use the inDrop1 batch, which includes"
    "16382 cells × 18771 genes.",
    image="openproblems-python-extras",
)
def pancreas(test=False):
    adata = load_pancreas(test=test)
    adata = adata[adata.obs["tech"].isin(["inDrop1"])]
    adata = utils.split_data(adata)
    return adata
