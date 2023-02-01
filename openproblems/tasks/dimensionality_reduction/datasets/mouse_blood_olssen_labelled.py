from ....data.mouse_blood_olssen_labelled import load_olsson_2016_mouse_blood
from ....tools.decorators import dataset
from ....tools.normalize import log_cpm


@dataset(
    "Mouse myeloid lineage differentiation",
    data_url=load_olsson_2016_mouse_blood.metadata["data_url"],
    data_reference=load_olsson_2016_mouse_blood.metadata["data_reference"],
    dataset_summary=(
        "Myeloid lineage differentiation from mouse blood. "
        "Sequenced by SMARTseq in 2016 by Olsson et al. "
        "660 cells x 112815 features with 4 cell type labels"
    ),
)
def olsson_2016_mouse_blood(test=False):
    adata = load_olsson_2016_mouse_blood(test=test)
    adata.uns["n_genes"] = adata.shape[1]
    return log_cpm(adata)
