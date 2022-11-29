from .....tools.decorators import metric

"""
The kBET algorithm (v.0.99.6, release 4c9dafa) determines whether the label composition
of a k nearest neighborhood of a cell is similar to the expected (global) label
composition (Buettner et al., Nat Meth 2019). The test is repeated for a random subset
of cells, and the results are summarized as a rejection rate over all tested
neighborhoods. Thus, kBET works on a kNN graph.

We compute kNN graphs where k = 50 for joint embeddings and corrected feature outputs
via Scanpy preprocessing steps. To test for technical effects and to account for
cell-type frequency shifts across datasets, we applied kBET
separately on the batch variable for each cell identity label. Using the kBET defaults,
a k equal to the median of the number of cells per batch within each label is used for
this computation. Additionally, we set the minimum and maximum thresholds of k to 10 and
100, respectively. As kNN graphs that have been subset by cell identity labels may no
longer be connected, we compute kBET per connected component. If >25% of cells were
assigned to connected components too small for kBET computation (smaller than k × 3),
we assigned a kBET score of 1 to denote poor batch removal. Subsequently, kBET scores
for each label were averaged and subtracted from 1 to give a final kBET score.

In Open Problems we do not run kBET on graph outputs to avoid computation-intensive
diffusion processes being run.
"""


@metric(
    metric_name="kBET",
    paper_reference="luecken2022benchmarking",
    maximize=True,
    image="openproblems-r-extras",
)
def kBET(adata):
    from scib.metrics import kBET

    return kBET(
        adata,
        batch_key="batch",
        label_key="labels",
        type_="embed",
        embed="X_emb",
        scaled=True,
        verbose=False,
    )
