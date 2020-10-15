import numpy as np

from ....tools.decorators import metric


@metric(
    metric_name="Template metric",
    maximize=True,
    # image="openproblems-template-image" # only if required
)
def template_metric(adata):
    # TODO: update
    return np.mean(adata.obs["template_variable"] - adata.obs["template_output"])
