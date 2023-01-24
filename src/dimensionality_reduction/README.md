
# Dimensionality reduction for visualization

Reduction of high-dimensional datasets to 2D for visualization &
interpretation

## Task description

Dimensionality reduction is one of the key challenges in single-cell
data representation. Routine single-cell RNA sequencing (scRNA-seq)
experiments measure cells in roughly 20,000-30,000 dimensions (i.e.,
features - mostly gene transcripts but also other functional elements
encoded in mRNA such as lncRNAs). Since its inception,scRNA-seq
experiments have been growing in terms of the number of cells measured.
Originally, cutting-edge SmartSeq experiments would yield a few hundred
cells, at best. Now, it is not uncommon to see experiments that yield
over [100,000 cells](https://www.nature.com/articles/s41586-018-0590-4)
or even [\> 1 million cells](https://doi.org/10.1126/science.aba7721).

Each *feature* in a dataset functions as a single dimension. While each
of the \~30,000 dimensions measured in each cell contribute to an
underlying data structure, the overall structure of the data is
challenging to display in few dimensions due to data sparsity and the
[*“curse of
dimensionality”*](https://en.wikipedia.org/wiki/Curse_of_dimensionality)
(distances in high dimensional data don’t distinguish data points well).
Thus, we need to find a way to [dimensionally
reduce](https://en.wikipedia.org/wiki/Dimensionality_reduction) the data
for visualization and interpretation.

## More information

- [Benchmarking
  results](https://openproblems-experimental.netlify.app/results_v2/dimensionality_reduction/)
- [Documentation](https://openproblems-experimental.netlify.app/results_v2/dimensionality_reduction/documentation)
