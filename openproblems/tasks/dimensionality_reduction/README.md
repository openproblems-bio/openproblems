# Dimensionality reduction

## The task

Dimensionality reduction is one of the key challenges in single-cell data representation. Routine single-cell RNA sequencing (scRNA-seq) experiments measure cells in roughly 20,000-30,000 dimensions (i.e., features - mostly gene transcripts but also other functional elements encoded in mRNA such as lncRNAs.) Since its inception, scRNA-seq experiments have been growing in terms of the number of cells measured. Originally, cutting-edge SmartSeq experiments would yield <a href="#">a few hundred cells</a>, at best. Now, it is not uncommon to see experiments that yield over <a href="https://www.nature.com/articles/s41586-018-0590-4">100,000 cells</a> or even <a href="https://www.10xgenomics.com/blog/our-13-million-single-cell-dataset-is-ready-to-download">> 1 million cells</a>.

Each *feature* in a dataset functions as a single dimension. While each of the ~30,000 dimensions measured in each cell (not accounting for roughly 75-90% data dropout per cell, another issue entirely), likely contribute to some sort of data structure, the overall structure of the data is diluted due to the <a href = "https://en.wikipedia.org/wiki/Curse_of_dimensionality">*"curse of dimensionality"*</a>. In short, it's difficult to visualize the contribution of each individual gene in a way that makes sense to the human eye, i.e., two or three dimensions (at most). Thus, we need to find a way to <a href = "https://en.wikipedia.org/wiki/Dimensionality_reduction">dimensionally reduce</a> the data for visualization and interpretation.

## The metrics

### Root mean square error
---

$$
    RMSE = \sqrt{ \sum_{i=1}^{n} \frac{(\hat{y}_i - y_i)^2}{n} }
$$

Where $y_i$ is the sum of pairwise euclidian distances between each value embedded in low-dimensional space and $\hat{y_i}$ is the sum of pairwise euclidian distances between each value in the original, high-dimensional space. The goal, in terms of preservation of this space is to minimize the difference between these terms. Finding the root-mean of the square of all differences (Root mean square error or $RMSE$ is a simple way to represent this as a scalar, which can then be used to compare to other methods.

<a href = "http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964a.pdf">Kruskel's stress</a> uses the RMSE, more or less in the now commonly-used MDS (multi-dimensional scaling). We can calculate and plot Kruskel's stress to get an idea where the majority of distortion of the ***topography*** of the data in high-dimensional space.

### Trustworthiness
---

_Adapted from the [sklearn documentation.](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.trustworthiness.html)_

Trustworthiness expresses to what extent the local structure in an embedding is retained. The trustworthiness is within [0, 1]. It is defined as

$$
    T(k) = 1 - \frac{2}{nk (2n - 3k - 1)} \sum^n_{i=1}
    \sum_{j \in \mathcal{N}_{i}^{k}} \max(0, (r(i, j) - k)) 
$$

where for each sample i, $\mathcal{N}_{i}^{k}$ are its k nearest neighbors in the output space, and every sample j is its $r(i, j)$-th nearest neighbor in the input space. In other words, any unexpected nearest neighbors in the output space are penalised in proportion to their rank in the input space.

References:  
    * "Neighborhood Preservation in Nonlinear Projection Methods: An Experimental Study" J. Venna, S. Kaski  
    * "Learning a Parametric Embedding by Preserving Local Structure" L.J.P. van der Maaten  

## The results

![Stress-plot RMASE CITE-seq](https://i.imgur.com/ulyAF9j.png)


Above is a "complex heatmap", which aims to show the regions that contribute the most stress. You can see that while a majority of the stress comes from the left side of the plot (as shown by the top of the complex heat map), the center of that left set of clusters does not contribute much to the stress, leading us to believe that by the measure of RMSE, the topology is relatively well-preserved. The stress mostly comes from the clusters at the top and bottom of that group of clusters spread across the second PC.

![RMSE PCA CITE-seq](https://i.imgur.com/H5rvIf6.png)

We performed principle component analysis, obtaining the first 50 components. We can then calculate the relative stress using the RMSE for each, in comparison to the original data, $y$. As one might suspect, the more components used, the lower the amount of distortion of the original data.

![RMSE Comparison](https://i.imgur.com/EoI72TI.png)

We can make this comparison across multiple dimensionality reduction methods. We can see that *t*-SNE seems to distort the data the least, in terms of pairwise euclidian distances. This does not necessarily mean the data is best represented by *t*-SNE, however. There are multiple means of measuring the "goodness" of a dimensionality reduction; RMSE is simply one of them.

## API

**Methods** should assign dimensionally-reduced 2D embedding coordinates to `adata.obsm['X_emb']` as well as the method of choice, i.e., `adata.obsm['X_pca']` or `adata.obsm['X_tsne']`.

**Metrics** should calculate the quality or "goodness of fit" of a dimensional reduction **method**.

## Pre-processing

There is a `preprocessing.py` function which contains a set of standard pre-processing functions. These perform steps such as transforming the raw data, selecting features and summarising dimensionality reduction. Where possible each **method** should first call one of these functions and use the processed `adata.X` or `adata.obsm['X_input']` slot as the input to the method. Variants of methods can be created by applying different pre-processing prior to the method itself (see `phate.py` for an example).

1. Raimundo, F., Vallot, C. & Vert, J. Tuning parameters of dimensionality reduction methods for single-cell RNA-seq analysis. Genome Biol 21, 212 (2020). https://doi.org/10.1186/s13059-020-02128-7
