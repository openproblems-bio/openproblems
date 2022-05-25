from .....tools.decorators import metric


@metric(
    metric_name="kBET",
    maximize=True,
    image="openproblems-r-extras",  # only if required
)
def kBET(adata):
    from scib.metrics import kBET

    import numpy as np

    kbet_score = kBET(
        adata,
        batch_key="batch",
        label_key="labels",
        type_="embed",
        embed="X_emb",
        scaled=True,
        verbose=False,
    )
