from .....tools.decorators import metric

"""
The cell-cycle conservation score evaluates how well the cell-cycle effect can be
captured before and after integration. We computed cell-cycle scores using Scanpy’s
score_cell_cycle function with a reference gene set from Tirosh et al for the
respective cell-cycle phases. We used the same set of cell-cycle genes for mouse and
human data (using capitalization to convert between the gene symbols). We then computed
the variance contribution of the resulting S and G2/M phase scores using principal
component regression (Principal component regression), which was performed for each
batch separately. The differences in variance before, Varbefore, and after, Varafter,
integration were aggregated into a final score between 0 and 1, using the equation:
CCconservation=1−|Varafter−Varbefore|/Varbefore.

In this equation, values close to 0 indicate lower conservation and 1 indicates complete
conservation of the variance explained by cell cycle. In other words, the variance
remains unchanged within each batch for complete conservation, while any deviation from
the preintegration variance contribution reduces the score."""


@metric(
    metric_name="Cell Cycle Score",
    maximize=True,
    image="openproblems-python-batch-integration",
)
def cc_score(adata, test=False):
    from ._utils import _get_split
    from scib.metrics import cell_cycle

    try:
        cc = cell_cycle(
            *_get_split(adata), "batch", embed="X_emb", organism=adata.uns["organism"]
        )

    except ValueError:
        cc = 0
    return cc
