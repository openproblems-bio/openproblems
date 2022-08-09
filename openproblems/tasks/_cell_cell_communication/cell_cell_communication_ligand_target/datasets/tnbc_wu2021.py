from .....data.tnbc_wu2021 import load_tnbc_data
from .....tools.decorators import dataset


@dataset(
    "Triple negative breast cancer atlas",
    data_url=load_tnbc_data.metadata["data_url"],
    data_reference=load_tnbc_data.metadata["data_reference"],
    dataset_summary="A single-cell atlas of human breast cancers with inferred "
    "cytokine activities as assumed true cell-cell communication. Cytokine "
    "activities were estimated by fitting a multivariate linear model with "
    "cytokine-focused signatures (see Dimitrov et al., 2022).",
    image="openproblems-r-extras",
)
def tnbc_data(test=False):
    return load_tnbc_data(test=test)
