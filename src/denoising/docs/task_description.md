## The task

Single-cell RNA-Seq protocols only detect a fraction of the mRNA molecules present
in each cell. As a result, the measurements (UMI counts) observed for each gene and each
cell are associated with generally high levels of technical noise ([Grün et al.,
2014](https://www.nature.com/articles/nmeth.2930)). Denoising describes the task of
estimating the true expression level of each gene in each cell. In the single-cell
literature, this task is also referred to as *imputation*, a term which is typically
used for missing data problems in statistics. Similar to the use of the terms "dropout",
"missing data", and "technical zeros", this terminology can create confusion about the
underlying measurement process ([Sarkar and Stephens,
2020](https://www.biorxiv.org/content/10.1101/2020.04.07.030007v2)).

A key challenge in evaluating denoising methods is the general lack of a ground truth. A
recent benchmark study ([Hou et al.,
2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x))
relied on flow-sorted datasets, mixture control experiments ([Tian et al.,
2019](https://www.nature.com/articles/s41592-019-0425-8)), and comparisons with bulk
RNA-Seq data. Since each of these approaches suffers from specific limitations, it is
difficult to combine these different approaches into a single quantitative measure of
denoising accuracy. Here, we instead rely on an approach termed molecular
cross-validation (MCV), which was specifically developed to quantify denoising accuracy
in the absence of a ground truth ([Batson et al.,
2019](https://www.biorxiv.org/content/10.1101/786269v1)). In MCV, the observed molecules
in a given scRNA-Seq dataset are first partitioned between a *training* and a *test*
dataset. Next, a denoising method is applied to the training dataset. Finally, denoising
accuracy is measured by comparing the result to the test dataset. The authors show that
both in theory and in practice, the measured denoising accuracy is representative of the
accuracy that would be obtained on a ground truth dataset.

## The metrics

Metrics for data denoising aim to assess denoising accuracy by comparing the denoised
*training* set to the randomly sampled *test* set.

* **MSE**: The mean squared error between the denoised counts of the training dataset
  and the true counts of the test dataset after reweighting by the train/test ratio.
* **Poisson**: The Poisson log likelihood of observing the true counts of the test
  dataset given the distribution given in the denoised dataset.

## API

Datasets should contain the raw UMI counts in `adata.X`, subsampled to training
(`adata.obsm["train"]`) and testing (`adata.obsm["test"]`) datasets using
`openproblems.tasks.denoising.datasets.utils.split_data`.

The task-specific data loader functions should split the provided raw UMI counts into a
training and a test dataset, as described by [Batson et al.,
2019](https://www.biorxiv.org/content/10.1101/786269v1). The training dataset should be
stored in `adata.obsm['train']`, and the test dataset should be stored in
`adata.obsm['test']`. Methods should store the denoising result in
`adata.obsm['denoised']`. Methods should not edit `adata.obsm["train"]` or
`adata.obsm["test"]`.