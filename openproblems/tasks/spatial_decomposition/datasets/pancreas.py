from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset
from ._utils import generate_synthetic_dataset


@dataset("Pancreas (average)")
def pancreas_average(test=False):
    adata = load_pancreas(test=test)
    adata.obs["label"] = adata.obs["celltype"]

    adata_spatial = generate_synthetic_dataset(adata, sim_type="avg")
    return adata_spatial


@dataset("Pancreas (cell)")
def pancreas_cell(test=False):
    adata = load_pancreas(test=test)
    adata.obs["label"] = adata.obs["celltype"]

    adata_spatial = generate_synthetic_dataset(adata, sim_type="cell")
    return adata_spatial
