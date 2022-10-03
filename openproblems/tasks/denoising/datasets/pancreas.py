from .....data.pancreas import load_pancreas
from .....tools.decorators import dataset

import scanpy as sc


@dataset(
    dataset_name="Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    data_reference=load_pancreas.metadata["data_reference"],
    dataset_summary="Human pancreatic islet scRNA-seq data from 6 datasets "
    "across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, "
    "and SMARTER-seq).",
    image="openproblems",
)
def pancreas(test=False, integer_only=True, techkeeps = [-3]):
    adata = load_pancreas(test)
    adata = utils.split_data(adata)
    return adata
