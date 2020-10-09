from ....tools.decorators import dataset
from ....data.scicar_chrom import load_scicar_mouse_kidney


@dataset("sciCAR Mouse Kidney with cell clusters")
def scicar_mouse_kidney(test=False):
    return load_scicar_mouse_kidney(test=test)
