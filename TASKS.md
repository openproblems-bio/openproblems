# Benchmarking Task Descriptions

Table of Contents
  * [Predicting gene expression from chromatin accessibility](#predicting-gene-expression-from-chromatin-accessibility)
  * [Data denoising and imputation](#data-denoising-and-imputation)
  * [Dimensionality reduction for visualization](#dimensionality-reduction-for-visualization)
  * [Cell type label projection from a reference atlas](#cell-type-label-projection-from-a-reference-atlas)
  * [Multimodal data integration (e.g. combining CITE-seq with RNA-seq)](#multimodal-data-integration--eg-combining-cite-seq-with-rna-seq-)
  * [Differential abundance from experimental perturbations](#differential-abundance-from-experimental-perturbations)
  * [Data integration and batch normalization](#data-integration-and-batch-normalization)
  * [Further tasks](#further-tasks)


### Predicting gene expression from chromatin accessibility

Several recently published techniques allow for simultaneous measurement of chromatin accessibility and gene expression in the single-cell level. For example, [SHARE-seq](https://pubmed.ncbi.nlm.nih.gov/33098772/) uses droplet-based sequencing to measure both modalities, [sci-CAR](https://doi.org/10.1126/science.aau0730) uses plate-based method to sequence both information. One challenge to analyze this data is to link the chromatin accessibility to gene expression in single-cells, which is useful to understand the single-cell trajectory and regulatory effect. 

Here, we compare methods for calculating gene scores (regulatory effect or chromatin potential) from the peak matrix of single-cell ATAC-seq modality, and evaluate their gene expression matrix prediction performance of single-cell RNA-seq by measuring consistency between the gene scores and gene expression across cells. 

Two representative software [ArchR](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html) and [MAESTRO](https://github.com/liulab-dfci/MAESTRO/blob/d58fc186a4329febde3b0d5b213c2d0edf8de44c/example/Gene_activity_modelling/Gene_activity_modelling.md), each has different way of assigning peak weights upon genes include a set of methods, but also share some common property. A typical algorithm includes the following common steps:

1) Given a peak matrix (cells x peaks), calculate the peaks weight matrix based on distances between gene transcription start site (TSS) and the peaks, this might involve reweighting potential Pol II elongation effect and masking of closest gene effect. 
2) Dot product between the peak matrix and peak weight matrix to get the gene score matrix (cells x genes). 
3) Evaluation of the regulatory effect in the single-cell levels by calculating the following metrics across genes for each single-cell, alternatively, the gene expression matrix and gene score matrix can be imputed by aggregating within nearest neighboring cells or cells within a cluster to get the pseudo-bulk gene expression matrix and pseudo-bulk gene score matrix (bulk/cluster x genes), and evaluate the metrics across genes for each bulk/cluster.

### Data denoising and imputation

Despite the proposed application of these protocols for describing per-cell information within a tissue or cell fraction of interest, single-cell RNA-Seq protocols only detect a small fraction of the mRNA molecules present in each cell. As a result, the measurements (UMI counts) observed for each gene and each cell are associated with generally high levels of technical noise ([Grün et al., 2014](https://www.nature.com/articles/nmeth.2930)). Denoising describes the task of estimating the true expression level of each gene in each cell. In the single-cell literature, this task is also referred to as imputation, a term which is typically used for missing data problems in statistics. Similar to the use of the terms "dropout", "missing data", and "technical zeros", this terminology can create confusion about the underlying measurement process ([Sarkar and Stephens, 2020](https://www.biorxiv.org/content/10.1101/2020.04.07.030007v2)).


A key challenge in evaluating denoising methods is the general lack of a ground truth. A recent benchmark study ([Hou et al., 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x)) relied on flow-sorted datasets, mixture control experiments ([Tian et al., 2019](https://www.nature.com/articles/s41592-019-0425-8)), and comparisons with bulk RNA-Seq data. Since each of these approaches suffers from specific limitations, it is difficult to combine these different approaches into a single quantitative measure of denoising accuracy. Here, we instead rely on an approach termed molecular cross-validation (MCV), which was specifically developed to quantify denoising accuracy in the absence of a ground truth ([Batson et al., 2019](https://www.biorxiv.org/content/10.1101/786269v1)). In MCV, the observed molecules in a given scRNA-Seq dataset are first partitioned between a training and a test dataset. Next, a denoising method is applied to the training dataset. Finally, denoising accuracy is measured by comparing the result to the test dataset. The authors show that both in theory and in practice, the measured denoising accuracy is representative of the accuracy that would be obtained on a ground truth dataset.
 
#### The metrics

Metrics for evaluating denoising of single cell datasets aim to characterize how well a denoised training subset compares to a held out testing subset, or otherwise to a secondary measurement (RNA-seq, single-cell atlas datasets, purified cell fractions). 

- MSE:
    - Let f(i)∈F be the denoised scRNA-seq measurement of cell i in the training set. Let g(j)∈G be the scRNA-seq measurement of cell j in the test set. MSE calculates the squared error difference between the observed cells in F and those in G. 
- Poisson Negative Log Likelihood Loss (Poisson NLLL):
    - Let f(i)∈F be the denoised scRNA-seq measurement of cell i in the training set. Let g(j)∈G be the scRNA-seq measurement of cell j in the test set. Poisson NLLL calculates the loss metric associated with likelihoods that F and G are sampled from the same Poisson distribution.

### Dimensionality reduction for visualization

Dimensional reduction is one of the key challenges in single-cell data representation. Routine single-cell RNA sequencing (scRNA-seq) experiments measure cells in roughly 20,000-30,000 dimensions (i.e., features - mostly gene transcripts but also other functional elements encoded in mRNA such as lncRNAs.) Since its inception, scRNA-seq experiments have been growing in terms of the number of cells measured. Originally, cutting-edge SmartSeq experiments would yield <a href="#">a few hundred cells</a>, at best. Now, it is not uncommon to see experiments that yield over <a href="https://www.nature.com/articles/s41586-018-0590-4">100,000 cells</a> or even <a href="https://www.10xgenomics.com/blog/our-13-million-single-cell-dataset-is-ready-to-download">> 1 million cells</a>.

Each *feature* in a dataset functions as a single dimension. While each of the ~30,000 dimensions measured in each cell (not accounting for roughly 75-90% data dropout per cell, another issue entirely), likely contribute to some sort of data structure, the overall structure of the data is diluted due to the <a href = "https://en.wikipedia.org/wiki/Curse_of_dimensionality">*"curse of dimensionality"*</a>. In short, it's difficult to visualize the contribution of each individual gene in a way that makes sense to the human eye, i.e., two or three dimensions (at most). Thus, we need to find a way to <a href = "https://en.wikipedia.org/wiki/Dimensionality_reduction">dimensionally reduce</a> the data for visualization and interpretation.

### Cell type label projection from a reference atlas

Label projection refers to the automatic identification of cell identity labels in a test or query dataset based on a reference dataset (or datasets) that typically contains curated, manually labeled cells. This process enables the rapid annotation of new datasets, which is becoming increasingly important as data generation becomes easier and larger in scale.

Label projection methods range from logistic regression which ignores batch information, to methods that perform projection based on a batch-integrated embedding. For a review, see [Abdelaal et al. (2019)](https://doi.org/10.1186/s13059-019-1795-z).

### Multimodal data integration (e.g. combining CITE-seq with RNA-seq)

Multimodal data integration refers to the task of combining together two datasets of different modalities of measurements (e.g., single-cell RNA sequencing and single-cell ATAC sequencing) on different observations of the same biological system. Integrating such measurements allows us to analyze the interaction between the different modalities, without requiring an explicitly joint measurement like [sci-CAR](https://doi.org/10.1126/science.aau0730) or [CITE-seq](https://doi.org/10.1038/nmeth.4380).

Here, we use such joint measurement modalities as a gold standard for multimodal data alignment; since we know the joint measurements are obtained from the same cells, a method which successfully realigns cells from one modality to their paired measurements from another modality (without prior knowledge of the shared cell identity) is known to be performant on the given dataset.

### Differential abundance from experimental perturbations 

An increasing number of single-cell studies are designed to compare reference single-cell atlases to phenotypes measured in disease, development, aging or upon experimental treatments. Perturbed cell states can be detected as shifts in abundance of cells on the transcriptional state space. 
Our goal is to evaluate the accuracy and sensitivity of approaches to quantify Differential Abundance (DA).

There are several challenges associated with this analysis, that will likely require different subtasks:

- Testing on cell clusters or populations defined _a priori_ can reduce the sensitivity to detect shifts in small cell subpopulations and when cellular states form a continuum
- Accounting for compositional biases: with relatively small numbers of cells profiled per experimental sample, the depletion of a single cell population leads to a relative increase in frequency of other populations. Models that test abundance changes of different populations independently may falsely detect these enrichment as true effects (see [Buettner et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.14.422688v2.full.pdf) for a detailed description of this bias). 
- Replicate measurements from the same biological condition can show substantial variance in cell abundances
- Technical and confounding sources of variation can introduce differences in cell numbers between biological samples. This is especially true for comparative studies on human tissues, which we expect will scale up in size and experimental design complexity in the near future.

**Datasets:** An additional challenge is the lack of a ground-truth. To overcome this problem, we will simulate condition labels on true single-cell datasets based by generating a smooth underlying condition probability. This will require defining minima and maxima of condition probability and interpolating smoothly between these points, with an approach such as [Inverse Distance Weighting](https://en.wikipedia.org/wiki/Inverse_distance_weighting).

**Methods:** This task is commonly approached by testing for differences in composition of predefined and discrete cell type clusters ([Haber et al. 2017](https://www.nature.com/articles/nature24489), [scDC](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3211-9), [scCODA](https://www.biorxiv.org/content/10.1101/2020.12.14.422688v2.full.pdf)), while methods quantify differences over KNN graphs ([MELD](https://www.nature.com/articles/s41587-020-00803-5), [Milo](https://www.biorxiv.org/content/10.1101/2020.11.23.393769v1.full), [DAseq](https://www.biorxiv.org/content/10.1101/711929v3)) or differentiation trajectories ([_condiments_](https://www.biorxiv.org/content/10.1101/2021.03.09.433671v1)).

**Metrics:** during the jamboree we will define which subtasks to prioritize. In general, we aim to include metrics that assess both the accuracy of the estimate of the effect size of differential abundance (such as the Mean Squared Error between the precited effect size for DA and the true condition probability) as well as the ability to detect subpopulations of cells that are similar in gene expression space and similarly enriched or depleted by an experimental perturbation (this is practically  useful to enable downstream investigation of condition-specific cell states).

### Data integration and batch normalization

Batch (or data) integration methods integrate datasets across batches that arise from various biological (e.g., tissue, location, individual, species) and technical (e.g., ambient RNA, lab, protocol) sources. The goal of a batch integration method is to remove unwanted batch effects in the data, while retaining biologically-meaningful variation that can help us to detect cell identities, fit cellular trajectories, or understand patterns of gene or pathway activity.

Methods that integrate batches typically have one of three different types of output: a corrected feature matrix, a joint embedding across batches, and/or an integrated cell-cell similarity graph (e.g., a kNN graph). In order to define a consistent input and output for each method and metric, we have divided the batch integration task into three subtasks. These subtasks are:

* Batch integration graphs
* Batch integration embeddings
* Batch integrated feature matrices

These subtasks collate methods that have the same data output type and metrics that evaluate this output. As corrected feature matrices can be turned into embeddings, which in turn can be processed into integrated graphs, methods overlap between the tasks. All methods are added to the graph subtask and imported into other subtasks from there.

Metrics for this task can be divided into those that assess the removal of batch effects, and assessments of the conservation of biological variation. This can be a helpful distinction when devising new metrics. This task, including the subtask structure, was taken from a [benchmarking study of data integration methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2). This is a useful reference for more background reading on this task and the above concepts.

### Further tasks

Next to the already initiated tasks, there are many open problems in computational single-cell research that we have not yet captured. These include spatial transcriptomics tasks (e.g., inferring spatial graphs, detecting spatially variable genes, deconvolving spots), perturbation modeling (e.g., predicting perturbed transcriptome distributions), or inferring gene regulatory networks. These tasks would benefit from a well-defined task definition and an associated standardized benchmarking setup. If you are motivated, please don’t hesitate to find a group of kindred spirits to initiate one of these tasks!
