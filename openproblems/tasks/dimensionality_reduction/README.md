# Dimensionality reduction

## The task

Dimensionality reduction is one of the key challenges in single-cell data
representation. Routine single-cell RNA sequencing (scRNA-seq) experiments measure cells
in roughly 20,000-30,000 dimensions (i.e., features - mostly gene transcripts but also
other functional elements encoded in mRNA such as lncRNAs). Since its inception,
scRNA-seq experiments have been growing in terms of the number of cells measured.
Originally, cutting-edge SmartSeq experiments would yield a few hundred cells, at best.
Now, it is not uncommon to see experiments that yield over [100,000 cells]
(<https://www.nature.com/articles/s41586-018-0590-4>) or even [> 1 million
cells.](https://www.10xgenomics.com/blog/our-13-million-single-cell-dataset-is-ready-to-download)

Each *feature* in a dataset functions as a single dimension. While each of the ~30,000
dimensions measured in each cell (not accounting for roughly 75-90% data dropout per
cell, another issue entirely), likely contribute to some sort of data structure, the
overall structure of the data is diluted due to the [*"curse of
dimensionality"*](https://en.wikipedia.org/wiki/Curse_of_dimensionality). In short, it's
difficult to visualize the contribution of each individual gene in a way that makes
sense to the human eye, i.e., two or three dimensions (at most). Thus, we need to find a
way to [dimensionally reduce](https://en.wikipedia.org/wiki/Dimensionality_reduction)
the data for visualization and interpretation.

## The methods

### Principal component analysis (PCA)

[Reference](https://www.tandfonline.com/doi/abs/10.1080/14786440109462720)

*Adapted from the [sklearn
documentation](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)*.

Linear dimensionality reduction using Singular Value Decomposition of the data to
project it to a lower dimensional space. The input data is centered but not scaled for
each feature before applying the SVD.

It uses the `scipy.sparse.linalg` ARPACK implementation of the truncated SVD as provided
by [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.pca.html).

### t-distributed stochastic neighbor embedding (t-SNE)

[Reference](https://arxiv.org/abs/1802.03426)

*Adapted from the [sklearn
documentation](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)*.

t-SNE is a tool to visualize high-dimensional data. It converts similarities between
data points to joint probabilities and tries to minimize the Kullback-Leibler divergence
between the joint probabilities of the low-dimensional embedding and the
high-dimensional data. t-SNE has a cost function that is not convex, i.e. with different
initializations we can get different results.

It is highly recommended to use another dimensionality reduction method to reduce the
number of dimensions to a reasonable amount (e.g. 50) if the number of features is very
high. This will suppress some noise and speed up the computation of pairwise distances
between samples. The implemented version first applies PCA with 50 dimensions before
calling the function provided by
[scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.tsne.html).

### Uniform manifold approximation and projection (UMAP)

[Reference](https://arxiv.org/abs/1802.03426)

*Adapted from the [umap-learn
documentation](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html)*.

UMAP is an algorithm for dimension reduction based on manifold learning techniques and
ideas from topological data analysis. The first phase consists of constructing a fuzzy
topological representation, based on nearest neighbours. The second phase is simply
optimizing the low dimensional representation to have as close a fuzzy topological
representation as possible to the full-dimensional data as measured by cross entropy.

The implemented version first applies PCA with 50 dimensions and calculates a
nearest-neighbour graph before calling the modified implementation in
[scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html).

### densMAP

[Reference](https://www.nature.com/articles/s41587-020-00801-7)

A modification of UMAP that adds an extra cost term in order to preserve information
about the relative local density of the data.

The implemented version first applies PCA with 50 dimensions before calling the function
from [umap-learn](https://umap-learn.readthedocs.io/en/latest/densmap_demo.html).

**Variants:**

* The (logCPM-normalized, 1000 HVG) expression matrix
* 50 principal components

### Potential of heat-diffusion for affinity-based transition embedding (PHATE)

[Reference](https://www.nature.com/articles/s41587-019-0336-3)

The five main steps of PHATE are:

1. Compute the pairwise distances from the data matrix.
2. Transform the distances to affinities to encode local information.
3. Learn global relationships via the diffusion process.
4. Encode the learned relationships using the potential distance.
5. Embed the potential distance information into low dimensions for visualization.

This implementation is from the [phate package](https://phate.readthedocs.io/en/stable/)

**Variants:**

* The square-root CPM transformed expression matrix
* 50 principal components of the logCPM-normalised, 1000 HVG expression matrix

### ivis

[Reference](https://www.nature.com/articles/s41598-019-45301-0)

[ivis](https://bering-ivis.readthedocs.io/en/latest/) is a machine learning library for
reducing dimensionality of very large datasets using Siamese Neural Networks.

### NeuralEE

[Reference](https://www.frontiersin.org/articles/10.3389/fgene.2020.00786/full)

A neural network implementation of elastic embedding implemented in the [NeuralEE
package](https://neuralee.readthedocs.io/en/latest/).

**Variants:**

* Scaled 500 HVGs from a logged expression matrix (no library size normalization)
* LogCPM-normalised, 1000 HVG expression matrix

### scvis

[Reference](https://www.nature.com/articles/s41467-018-04368-5#Sec10)

A neural network generative model that uses the t-SNE objective as a constraint
implemented in the [scvis package](https://bitbucket.org/jerry00/scvis-dev/).

## The metrics

### Root mean square error

$$
    RMSE = \sqrt{ \sum_{i=1}^{n} \frac{(\hat{y}_i - y_i)^2}{n} }
$$

Where $y_i$ is the sum of pairwise euclidean distances between each value embedded in
low-dimensional space and $\hat{y_i}$ is the sum of pairwise euclidean distances between
each value in the original, high-dimensional space. The goal, in terms of preservation
of this space is to minimize the difference between these terms. Finding the root-mean
of the square of all differences (Root mean square error or $RMSE$ is a simple way to
represent this as a scalar, which can then be used to compare to other methods.

[Kruskel's
stress](http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964a.pdf)
uses the RMSE, more or less in the now commonly-used MDS (multi-dimensional scaling). We
can calculate and plot Kruskel's stress to get an idea where the majority of distortion
of the ***topography*** of the data in high-dimensional space.

### Trustworthiness

*Adapted from the [sklearn
documentation.](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.trustworthiness.html)*

Trustworthiness expresses to what extent the local structure in an embedding is
retained. The trustworthiness is within [0, 1]. It is defined as

$$
    T(k) = 1 - \frac{2}{nk (2n - 3k - 1)} \sum^n_{i=1}
    \sum_{j \in \mathcal{N}_{i}^{k}} \max(0, (r(i, j) - k))
$$

where for each sample i, $\mathcal{N}_{i}^{k}$ are its k nearest neighbors in the output
space, and every sample j is its $r(i, j)$-th nearest neighbor in the input space. In
other words, any unexpected nearest neighbors in the output space are penalised in
proportion to their rank in the input space.

References:

* "Neighborhood Preservation in Nonlinear Projection Methods: An Experimental Study" J.
  Venna, S. Kaski
* "Learning a Parametric Embedding by Preserving Local Structure" L.J.P. van der Maaten

### Density preservation

[Reference](https://www.nature.com/articles/s41587-020-00801-7)

The local density preservation metric that is part of the cost function for densMAP.
Some parts of this are re-implemented as the are not exposed my the [umap-learn
package.](https://umap-learn.readthedocs.io)

### NN Ranking

[Reference](https://doi.org/10.1016/j.heliyon.2021.e06199)

A set of metrics from the [pyDRMetrics package](https://doi.org/10.17632/jbjd5fmggh.2).
The implementation uses a slightly modified version of the original source code rather
than the PyPI package that is now available.

* Continuity - Measures error of hard extrusions
* co-KNN size - Counts how many points are in both k-nearest neighbors before and after
  the DR
* co-KNN AUC - The area under the co-KNN curve
* Local continuity meta criterion - co-KNN size with baseline removal which favors
  locality more
* Local property metric - Summary of the local co-KNN
* Global property metric - Summary of the global co-KNN

## Example results

![Stress-plot RMASE CITE-seq](https://i.imgur.com/ulyAF9j.png)

Above is a "complex heatmap", which aims to show the regions that contribute the most
stress. You can see that while a majority of the stress comes from the left side of the
plot (as shown by the top of the complex heat map), the center of that left set of
clusters does not contribute much to the stress, leading us to believe that by the
measure of RMSE, the topology is relatively well-preserved. The stress mostly comes from
the clusters at the top and bottom of that group of clusters spread across the second
PC.

![RMSE PCA CITE-seq](https://i.imgur.com/H5rvIf6.png)

We performed principle component analysis, obtaining the first 50 components.
We can then calculate the relative stress using the RMSE for each, in comparison to the
original data, $y$. As one might suspect, the more components used, the lower the amount
of distortion of the original data.

![RMSE Comparison](https://i.imgur.com/EoI72TI.png)

We can make this comparison across multiple dimensionality reduction methods. We can see
that *t*-SNE seems to distort the data the least, in terms of pairwise euclidean
distances. This does not necessarily mean the data is best represented by *t*-SNE,
however. There are multiple means of measuring the "goodness" of a dimensionality
reduction; RMSE is simply one of them.

## API

**Datasets** should provide un-normalized raw counts in `adata.X`.

**Methods** should assign dimensionally-reduced 2D embedding coordinates to
`adata.obsm['X_emb']`.

**Metrics** should calculate the quality or "goodness of fit" of a dimensional reduction
**method**. If the un-normalized input counts matrix is required by the matrix it can be
accessed in `adata.layers["counts"]`.

## Pre-processing

Different methods can require different pre-processing of the data. Standard
pre-processing functions are available as part of the `tools` module. Where possible
each **method** should first call one of these functions and use the processed `adata.X`
slot as the input to the method. Raw counts are also stored in `adata.layers["counts"]`
by the standard pre-processing functions, if a method performs its own pre-processing it
should also do this for use by metrics. For most methods a standard pre-processing with
the `log_cpm_hvg()` function is used which normalizes the expression matrix to counts
per million (CPM), performs a log transformation and subsets the data to highly-variable
genes (HVGs) as selected by scanpy's `high_variable_genes(adata, n_top_genes=n_genes,
flavor="cell_ranger")` (1000 genes by default). Variants of methods can be created by
applying different pre-processing prior to the method itself (see `phate.py` for an
example).

1. Raimundo, F., Vallot, C. & Vert, J. Tuning parameters of dimensionality reduction
   methods for single-cell RNA-seq analysis. Genome Biol 21, 212 (2020).
   <https://doi.org/10.1186/s13059-020-02128-7>
