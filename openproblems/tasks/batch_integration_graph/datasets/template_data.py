import numpy as np
from ....data.template import load_template_data
from ....tools.decorators import dataset


@dataset("Template dataset")
def template_data(test=False):
    # TODO: update
    adata = load_template_data(test=test)
    # Apply any further processing
    adata.obs["template_variable"] = np.random.normal(0, 1, adata.shape[0])
    return adata
