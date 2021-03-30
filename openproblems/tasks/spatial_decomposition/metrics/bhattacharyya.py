from ....tools.decorators import metric

import numpy as np
import sklearn.metrics


@metric(metric_name="bhattacharyya coefficient")
