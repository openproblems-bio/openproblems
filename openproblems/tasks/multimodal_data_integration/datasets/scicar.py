from ....data.scicar import load_scicar_cell_lines, load_scicar_mouse_kidney


def scicar_cell_lines(test=False):
    return load_scicar_cell_lines(test=test)


def scicar_mouse_kidney(test=False):
    return load_scicar_mouse_kidney(test=test)
