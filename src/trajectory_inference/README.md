# Trajectory inference

## Task description
Trajectory inference (TI) is a computational analysis used in single-cell transcriptomics to determine the pattern of a dynamic process experienced by cells and then arrange cells based on their progression through the process. 
A trajectory is a graph where the nodes represent noteworthy cellular states, and each cell is predicted to be progressing along transitions between the different states (Figure 1A).
Main applications of TI are identifying branch points, end states, predicting the topology of the dynamic process, or identifying genes whose expression varies gradually along the topology (Figure 1B).

| ![](docs/images/trajectory_inference.png) | 
|:--:| 
| **Figure 1**: Trajectory inference for single-cell omics data. Image borrowed from [1]. **A**: During a dynamic process cells pass through several transitional states, characterized by different waves of transcriptional, morphological, epigenomic and/or surface marker changes [2]. TI methods provide an unbiased approach to identifying and correctly ordering different transitional stages. **B**: By overlaying gene expression levels on a dimensionality reduction, the milestones can be annotated to allow better interpretation of the cellular heterogeneity. |

A comparison of 45 TI methods on 110 real and 229 synthetic datasets found that the different methods are very complementary when comparing different types of input datasets, and that performance of a method can be highly variable even in multiple runs on the same input dataset [3]. 

A persisting issue amongst TI methods is the usage of a standard definition of the task and usage of well-defined input and output data structures in order to make results comparable between methods. This task uses more restrictive version of the data structures proposed by Saelens et al. [3], but updated to make use of the anndata file format as used in the rest of the openproblems project.

## Metrics

## API



## References
1. Robrecht Cannoodt. “Modelling single-cell dynamics with trajectories and gene regulatory networks“. Doctoral dissertation, Ghent University (2019). URL: [cannoodt.dev/files/phdthesis.pdf](https://cannoodt.dev/files/phdthesis.pdf).

2. Tariq Enver et al. “Stem Cell States, Fates, and the Rules of Attraction”. In: Cell Stem Cell 4.5 (May 8, 2009), pp. 387–397. ISSN: 1875-9777. DOI: [10.1016/j.stem.2009.04.011](https://doi.org/10.1016/j.stem.2009.04.011). pmid: 19427289.

3. Wouter Saelens, Cannoodt Robrecht et al. “A Comparison of Single-Cell Trajectory Inference Methods“. In: Nature Biotechnology 37 (May 2019). ISSN: 15461696. DOI: [10.1038/s41587-019-0071-9](https://doi.org/10.1038/s41587-019-0071-9).