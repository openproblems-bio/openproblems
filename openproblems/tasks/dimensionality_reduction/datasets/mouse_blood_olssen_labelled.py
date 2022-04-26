from ....data.mouse_blood_olssen_labelled import load_olsson_2016_mouse_blood
from ....tools.decorators import dataset


@dataset("Mouse blood. Olsson, et al. Nature. 2016")
def olsson_2016_mouse_blood(test=False):
    return load_olsson_2016_mouse_blood(test=test)
