# Batch integration

Batch (or data) integration methods integrate datasets across batches that arise from
various biological (e.g., tissue, location, individual, species) and technical (e.g.,
ambient RNA, lab, protocol) sources. The goal of a batch integration method is to remove
unwanted batch effects in the data, while retaining biologically-meaningful variation
that can help us to detect cell identities, fit cellular trajectories, or understand
patterns of gene or pathway activity.

Methods that integrate batches typically have one of three different types of output:
a corrected feature matrix, a joint embedding across batches, and/or an integrated
cell-cell similarity graph (e.g., a kNN graph). In order to define a consistent input
and output for each method and metric, we have divided the batch integration task into
three subtasks. These subtasks are:

* [Batch integration graphs](batch_integration_graph/),
* [Batch integration embeddings](batch_integration_embed/), and
* [Batch integrated feature matrices](batch_integration_feature/)

These subtasks collate methods that have the same data output type and metrics that
evaluate this output. As corrected feature matrices can be turned into embeddings, which
in turn can be processed into integrated graphs, methods overlap between the tasks. All
methods are added to the graph subtask and imported into other subtasks from there.
Information on the task API for datasets, methods, and metrics can be found in the
individual subtask pages.

Metrics for this task can be divided into those that assess the removal of batch
effects, and assessments of the conservation of biological variation. This can be a
helpful distinction when devising new metrics. This task, including the subtask
structure, was taken from a [benchmarking study of data integration
methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2). This is a useful
reference for more background reading on the task and the above concepts.
