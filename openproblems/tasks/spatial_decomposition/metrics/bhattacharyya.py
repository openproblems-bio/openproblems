from ....tools.decorators import metric

import sklearn.metrics
import numpy as np

@metric(metric_name="bhattacharyya coefficient")
