from ....tools.decorators import metric

import scIB

@metric(metric_name="Normalised mutual information", maximize=True)
def nmi(adata, method="arithmetic", nmi_dir=None):
    return scIB.metrics.nmi(adata, 'labels', 'clusters_post', method, nmi_dir)
