# Cell-cell Communication

## The task

The growing availability of single-cell data has sparked an increased
interest in the inference of cell-cell communication (CCC),
with an ever-growing number of computational tools developed for this purpose.

Different tools propose distinct preprocessing steps with diverse
scoring functions, that are challenging to compare and evaluate.
Furthermore, each tool typically comes with its own set of prior knowledge.
To harmonize these, [Dimitrov et
al, 2022](https://doi.org/10.1038/s41467-022-30755-0) recently developed the
[LIANA](https://github.com/saezlab/liana) framework, which was used
as a foundation for this task.

The challenges in evaluating the tools are further exacerbated by the
lack of aa gold standard to benchmark the performance of CCC methods. In an
attempt to address this, Dimitrov et al use alternative data modalities, including
the spatial proximity of cell types and inferred
downstream cytokine activities, to generate an inferred ground truth. However,
these modalities are only approximations of biological reality and come
with their own assumptions and limitations. In time, the inclusion of more
datasets with known ground truth interactions will become available, from
which the limitations and advantages of the different CCC methods will
be better understood.

This subtask evaluates methods that predict interations between ligand molecules and
target cell types.

## The metrics

Metrics for cell-cell communication aim to characterize how good are
the different scoring methods at prioritizing assumed truth predictions.

* **Odds ratio**: The odds ratio represents the ratio of true and false
positives within a set of prioritized interactions (top ranked hits) versus
the same ratio for the remainder of the interactions. Thus, in this
scenario odds ratios quantify the strength of association between the
ability of methods to prioritize interactions and those interactions
assigned to the positive class.

## API

### Datasets

Datasets should include cell type annotations in `adata.obs["label"]`, along with some
assumed truth in `adata.uns["ccc_target"]`. The assumed truth could be derived from
various proxies; we refer the reader to [Dimitrov et
al](https://doi.org/10.1038/s41467-022-30755-0) for more details.

`adata.uns["ccc_target"]` should be a Pandas DataFrame containing all of the following
columns:

* `response`: `int`, binary response variable indicating whether an interaction is
  assumed to have occurred
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `target`: `str`, name of target cell type in interaction

The datasets should also include a
[NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
in `adata.uns["target_organism"]` - used to convert the (typically human) prior
knowledge of the CCC methods to the corresponding gene homologs.
`adata.X` should contain the raw counts matrix.

TODO: describe ligand-receptor resource.

### Methods

Methods should predict interactions between cell types without using
`adata.uns["ccc_target"]`. Predicted interactions should be stored in
`adata.uns["ccc_pred"]` as a Pandas DataFrame containing all of the following columns:

* `score`: `float`, score between `-inf` to `+inf` giving a predicted strength of the
  inferred interaction
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `target`: `str`, name of target cell type in interaction

Methods should infer a score for each _intersecting interaction_ in the harmonized
prior-knowledge resource provided by LIANA. We define _intersecting interactions_ as
those for which the relevant genes are both present in the dataset and the resource.

The prior-knowledge resource is available via the
`cell_cell_communication.utils.ligand_receptor_resource` function, which returns a
DataFrame containing the columns `ligand_genesymbol` and `receptor_genesymbol`, which
correspond to the ligand and receptor genes, respectively. These may contain complexes
with subunits separated with `_`. Hence, **methods should be able to deal with
complex-containing interactions**.

### Metrics

Metrics should evaluate the concordance between `adata.uns["ccc_target"]` and
`adata.uns["ccc_pred"]` to evaluate the success of a method in predicting interactions.
