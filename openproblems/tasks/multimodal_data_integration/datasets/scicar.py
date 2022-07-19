from ....data.multimodal.scicar import load_scicar_cell_lines
from ....data.multimodal.scicar import load_scicar_mouse_kidney
from ....tools.decorators import dataset


@dataset(
    "sciCAR Cell Lines",
    data_url=load_scicar_cell_lines.metadata["data_url"],
    data_reference=load_scicar_cell_lines.metadata["data_reference"],
    dataset_summary="TODO",
)
def scicar_cell_lines(test=False):
    return load_scicar_cell_lines(test=test)


@dataset(
    "sciCAR Mouse Kidney",
    data_url=load_scicar_mouse_kidney.metadata["data_url"],
    dataset_summary="TODO",
)
def scicar_mouse_kidney(test=False):
    return load_scicar_mouse_kidney(test=test)
