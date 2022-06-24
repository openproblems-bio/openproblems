from .....tools.decorators import metric

"""
The graph connectivity metric assesses whether the kNN graph representation,
G, of the integrated data directly connects all cells with the same cell
identity label. For each cell identity label c, we created the subset kNN
graph G(Nc;Ec) to contain only cells from a given label. Using these subset
kNN graphs, we computed the graph connectivity score using the equation:

gc =1/|C| Σc∈C |LCC(G(Nc;Ec))|/|Nc|.

Here, C represents the set of cell identity labels, |LCC()| is the number
of nodes in the largest connected component of the graph, and |Nc| is the
number of nodes with cell identity c. The resultant score has a range
of (0;1], where 1 indicates that all cells with the same cell identity
are connected in the integrated kNN graph, and the lowest possible score
indicates a graph where no cell is connected. As this score is computed
on the kNN graph, it can be used to evaluate all integration outputs.
"""


@metric(
    metric_name="Graph connectivity",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def graph_connectivity(adata):
    import scib.metrics

    adata.obs["labels"] = adata.obs["labels"].astype("category")
    return scib.metrics.graph_connectivity(adata, "labels")
