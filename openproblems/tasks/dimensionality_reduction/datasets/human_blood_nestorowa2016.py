from ....data.human_blood_nestorowa2016 import load_human_blood_nestorowa2016
from ....tools.decorators import dataset


@dataset(
    "Mouse hematopoietic stem cell differentiation",
    data_url=load_human_blood_nestorowa2016.metadata["data_url"],
    data_reference=load_human_blood_nestorowa2016.metadata["data_reference"],
    dataset_summary="1.6k hematopoietic stem and progenitor cells from mouse bone "
    "marrow. Sequenced by Smart-seq2.",
)
def human_blood_nestorowa2016(test=False):
    return load_human_blood_nestorowa2016(test=test)
