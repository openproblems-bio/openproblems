# Matching modalities

## The task

Cellular function is regulated by the complex interplay of different types of biological
molecules (DNA, RNA, proteins, etc.), which determine the state of a cell. Several
recently described technologies allow for simultaneous measurement of different aspects
of cellular state. For example, [sci-CAR](https://openproblems.bio/bibliography#cao2018joint)
jointly profiles RNA expression and chromatin accessibility on the same cell and
[CITE-seq](https://openproblems.bio/bibliography#stoeckius2017simultaneous) measures
surface protein abundance and RNA expression from each cell. These technologies enable
us to better understand cellular function, however datasets are still rare and there are
tradeoffs that these measurements make for to profile multiple modalities.

Joint methods can be more expensive or lower throughput or more noisy than measuring a
single modality at a time. Therefore it is useful to develop methods that are capable
of integrating measurements of the same biological system but obtained using different
technologies on different cells.

Here the goal is to learn a latent space where cells profiled by different technologies in
different modalities are matched if they have the same state. We use jointly profiled
data as ground truth so that we can evaluate when the observations from the same cell
acquired using different modalities are similar. A perfect result has each of the paired
observations sharing the same coordinates in the latent space.

## The metrics

Metrics for matching modalities aim to characterize how well the aligned
datasets correspond to the ground truth.

* **kNN AUC**: Let $f(i) ∈ F$ be the modality 1 (e.g., scRNA-seq) measurement of cell $i$,
  and $g(i) ∈ G$ be the modality 2 (e.g., scATAC-seq) measurement of cell $i$. kNN-AUC
  calculates the average percentage overlap of neighborhoods of $f(i)$ in $F$ with
  neighborhoods of $g(i)$ in $G$. Higher is better.
* **MSE**: Mean squared error (MSE) is the average distance between each pair of matched
  observations of the same cell in the learned latent space. Lower is better.

## API

Datasets should include matched measurements from two modalities, which are contained in
`adata` and `adata.obsm["mode2"]`. The task is to align these two modalities as closely
as possible, without using the known bijection between the datasets.

Methods should create joint matrices `adata.obsm["aligned"]` and
`adata.obsm["mode2_aligned"]` which reside in a joint space.

Metrics should evaluate how well the cells which are known to be equivalent are aligned
in the joint space.
