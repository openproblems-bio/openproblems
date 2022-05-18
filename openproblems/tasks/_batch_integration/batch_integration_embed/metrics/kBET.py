from .....tools.decorators import metric


@metric(
    metric_name="kBET",
    maximize=True,
    image="openproblems-r-extras",  # only if required
)
def kBET(adata):
    from scib.metrics import kBET

    import numpy as np

    return 1 - np.nanmean(
        kBET(
            adata,
            batch_key="batch",
            label_key="label",
            embed="X_emb",
            subsample=0.5,
            heuristic=True,
            verbose=False,
        )["kBET"]
    )
