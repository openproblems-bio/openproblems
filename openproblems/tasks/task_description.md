# Dimensionality reduction

---

## The task

Dimensional reduction is one of the key challenges in single-cell data representation. Routine single-cell RNA sequencing (scRNA-seq) experiments measure cells in roughly 20,000-30,000 dimensions (i.e., features - mostly gene transcripts but also other functional elements encoded in mRNA such as lncRNAs.) Since its inception, scRNA-seq experiments have been growing in terms of the number of cells measured. Originally, cutting-edge SmartSeq experiments would yield <a href="#">a few hundred cells</a>, at best. Now, it is not uncommon to see experiments that yield over <a href="https://www.nature.com/articles/s41586-018-0590-4">100,000 cells</a> or even <a href="https://www.10xgenomics.com/blog/our-13-million-single-cell-dataset-is-ready-to-download">> 1 million cells</a>.

Each *feature* in a dataset functions as a single dimension. While each of the ~30,000 dimensions measured in each cell (not accounting for roughly 75-90% data dropout per cell, another issue entirely), likely contribute to some sort of data structure, the overall structure of the data is diluted due to the <a href = "https://en.wikipedia.org/wiki/Curse_of_dimensionality">*"curse of dimensionality"*</a>. In short, it's difficult to visualize the contribution of each individual gene in a way that makes sense to the human eye, i.e., two or three dimensions (at most). Thus, we need to find a way to <a href = "https://en.wikipedia.org/wiki/Dimensionality_reduction">dimensionally reduce</a> the data for visualization and interpretation.

## The metrics

### The root mean square error
---

$$
    RMSE = \sqrt{ \sum_{i=1}^{n} \frac{(\hat{y}_i - y_i)^2}{n} }
$$

Where $y_i$ is the sum of pairwise euclidian distances between each value embedded in low-dimensional space and $\hat{y_i}$ is the sum of pairwise euclidian distances between each value in the original, high-dimensional space. The goal, in terms of preservation of this space is to minimize the difference between these terms. Finding the root-mean of the square of all differences (Root mean square error or $RMSE$ is a simple way to represent this as a scalar, which can then be used to compare to other methods.

<a href = "http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964a.pdf">Kruskel's stress</a> uses the RMSE, more or less in the now commonly-used MDS (multi-dimensional scaling). We can calculate and plot Kruskel's stress to get an idea where the majority of distortion of the ***topography*** of the data in high-dimensional space.

## The results

![Stress-plot RMASE CITE-seq](https://i.imgur.com/ulyAF9j.png)


Above is a "complex heatmap", which aims to show the regions that contribute the most stress. You can see that while a majority of the stress comes from the left side of the plot (as shown by the top of the complex heat map), the center of that left set of clusters does not contribute much to the stress, leading us to believe that by the measure of RMSE, the topology is relatively well-preserved. The stress mostly comes from the clusters at the top and bottom of that group of clusters spread across the second PC.

![RMSE PCA CITE-seq](https://i.imgur.com/H5rvIf6.png)

We performed principle component analysis, obtaining the first 50 components. We can then calculate the relative stress using the RMSE for each, in comparison to the original data, $y$. As one might suspect, the more components used, the lower the amount of distortion of the original data.

![RMSE Comparison](https://i.imgur.com/EoI72TI.png)

We can make this comparison across multiple dimensionality reduction methods. We can see that *t*-SNE seems to distort the data the least, in terms of pairwise euclidian distances. This does not necessarily mean the data is best represented by *t*-SNE, however. There are multiple means of measuring the "goodness" of a dimensional reduction; RMSE is simply one of them.

# Data integration

---



# Regulatory effect prediction

Using chromatin accessibility to predict gene expression in single cells of multi-omics sequencing data for evaluation of regulatory effect

## The task
Several recently published techniques allow for simultaneous measurement of chromatin accessibility and gene expression in the single cell level. For example, SHARE-seq uses droplet-based sequencing to measure both modalities, sci-CAR uses plate-based method to sequence both information. One challenge to analyze this data is to link the chromatin accessibility to gene expression in single cells, which is useful to understand the single cell trajectory and regulatory effect.

Here, we compare methods for calculating gene scores (regulatory effect or chromatin potential) from the peak matrix of single cell ATAC-seq modality, and evaluate their gene expression matrix prediction performance of single cell RNA-seq by measuring consistency between the gene scores and gene expression across cells.

Two representative software [ArchR](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html) and [MAESTRO](https://github.com/liulab-dfci/MAESTRO/blob/d58fc186a4329febde3b0d5b213c2d0edf8de44c/example/Gene_activity_modelling/Gene_activity_modelling.md) include a set of methods, each has different way of assigning peak weights upon genes, but also share some common property. A typical algorithm includes the following common steps: 1) given a peak matrix (cells x peaks), calculate the peaks weight matrix based on distances between gene transcription start site (TSS) and the peaks, this might involve reweighting potential Pol II elongation effect and masking of closest gene effect. 2) Dot product between the peak matrix and peak weight matrix to get the gene score matrix (cells x genes). 3) evaluation of the regulatory effect in the single cell levels by calculating the following metrics across genes for each single cell, alternatively, the gene expression matrix and gene score matrix can be imputed by aggregating within nearest neighboring cells or cells within a cluster to get the pseudo-bulk gene expression matrix and pseudo-bulk gene score matrix (bulk/cluster x genes), and evaluate the metrics across genes for each bulk/cluster.

## The metrics
Metrics for regulatory effect prediction aim to characterize how well each gene score correctly predict the gene expression in the same cells or same clusters of cells.

Average Pearson/Spearman correlation: average correlation measurement between the gene scores and gene expression across cells or bulk/clusters of cells.


Mean squared error: average mean squared error measurement between the gene scores and gene expression across cells or bulk/clusters of cells.
