import numpy as np

from ....tools.decorators import metric


@metric(metric_name="Template metric", maximize=True)
def template_metric(adata):
    return np.mean(adata.obs["template_variable"] - adata.obs["template_output"])
