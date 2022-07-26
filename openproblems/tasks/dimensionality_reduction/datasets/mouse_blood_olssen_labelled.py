from ....data.mouse_blood_olssen_labelled import load_olsson_2016_mouse_blood
from ....tools.decorators import dataset


@dataset(
    "Mouse myeloid lineage differentiation",
    data_url=load_olsson_2016_mouse_blood.metadata["data_url"],
    data_reference=load_olsson_2016_mouse_blood.metadata["data_reference"],
    dataset_summary="Myeloid lineage differentiation from mouse blood. Sequenced by SMARTseq in 2016 by Olsson et al.",
)
def olsson_2016_mouse_blood(test=False):
    return load_olsson_2016_mouse_blood(test=test)
