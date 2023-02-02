from ....data.multimodal.scicar import load_scicar_cell_lines
from ....data.multimodal.scicar import load_scicar_mouse_kidney
from ....tools.decorators import dataset


@dataset(
    "sciCAR Cell Lines",
    data_url=load_scicar_cell_lines.metadata["data_url"],
    data_reference=load_scicar_cell_lines.metadata["data_reference"],
    dataset_summary=(
        "5k cells from a time-series of dexamethasone treatment sequenced "
        "by sci-CAR, a combinatorial indexing-based co-assay that jointly profiles "
        "chromatin accessibility and mRNA."
    ),
)
def scicar_cell_lines(test=False):
    return load_scicar_cell_lines(test=test)


@dataset(
    "sciCAR Mouse Kidney",
    data_url=load_scicar_mouse_kidney.metadata["data_url"],
    data_reference=load_scicar_cell_lines.metadata["data_reference"],
    dataset_summary=(
        "11k cells from adult mouse kidney sequenced "
        "by sci-CAR, a combinatorial indexing-based co-assay that jointly profiles "
        "chromatin accessibility and mRNA."
    ),
)
def scicar_mouse_kidney(test=False):
    return load_scicar_mouse_kidney(test=test)
