## The task

Several recently described technologies allow for simultaneous measurement of different aspects of cell state. For example, [sci-CAR](https://doi.org/10.1126/science.aau0730) jointly profiles RNA expression and chromatin accessibility on the same cell and [CITE-seq](https://doi.org/10.1038/nmeth.4380) measures surface protein abundance and RNA expression from each cell. However, these joint profiling methods have several tradeoffs compared to unimodal measurements.

Joint methods can be more expensive or lower throughput or more noisy than measuring a single modality at a time. Therefore it is useful to develop methods that are capable of integrating measurements of the same biological system but obtained using different technologies.

Here the goal is to learn a latent space where observations from the same cell acquired using different modalities. A perfect result has each of the paired observations sharing the same coordinates in the latent space.

## The metrics
Metrics for multimodal data integration aim to characterize how well the aligned datasets correspond to the ground truth.

* **kNN AUC**: Let $f(i) ∈ F$ be the scRNA-seq measurement of cell $i$, and $g(i) ∈ G$ be the scATAC- seq measurement of cell $i$. kNN-AUC calculates the average percentage overlap of neighborhoods of $f(i)$ in $F$ with neighborhoods of $g(i)$ in $G$. Higher is better.
* **MSE**: Mean squared error (MSE) is the average distance between each pair of matched observations of the same cell in the learned latent space. Lower is better.

## API

Datasets should include matched measurements from two modalities, which are contained in `adata` and `adata.obsm["mode2"]`. The task is to align these two modalities as closely as possible, without using the known bijection between the datasets.

Methods should create joint matrices `adata.obsm["aligned"]` and `adata.obsm["mode2_aligned"]` which reside in a joint space.

Metrics should evaluate how well the cells which are known to be equivalent are aligned in the joint space.
