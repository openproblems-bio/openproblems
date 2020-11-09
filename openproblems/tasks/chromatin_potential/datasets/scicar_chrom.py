from ....tools.decorators import dataset
from ....data.multimodal import scicar


@dataset("sciCAR Mouse Kidney with cell clusters")
def scicar_mouse_kidney(test=False):
    return scicar.load_scicar_mouse_kidney(test=test)
