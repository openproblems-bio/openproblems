from ....data.human_blood_nestorowa2016 import load_human_blood_nestorowa2016
from ....tools.decorators import dataset


@dataset(
    "Human blood (HSCs and differentiation thereof). Nestorowa, et al. Blood. 2016",
    data_url=load_human_blood_nestorowa2016.metadata["data_url"],
)
def human_blood_nestorowa2016(test=False):
    return load_human_blood_nestorowa2016(test=test)
