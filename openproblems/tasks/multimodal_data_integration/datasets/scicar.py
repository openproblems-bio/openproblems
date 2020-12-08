from ....data.multimodal.scicar import load_scicar_cell_lines
from ....data.multimodal.scicar import load_scicar_mouse_kidney
from ....tools.decorators import dataset


@dataset("sciCAR Cell Lines")
def scicar_cell_lines(test=False):
    return load_scicar_cell_lines(test=test)


@dataset("sciCAR Mouse Kidney")
def scicar_mouse_kidney(test=False):
    return load_scicar_mouse_kidney(test=test)
