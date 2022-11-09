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
lack of a gold standard to benchmark the performance of CCC methods. In an
attempt to address this, Dimitrov et al use alternative data modalities, including
the spatial proximity of cell types and 
downstream cytokine activities, to generate an inferred ground truth. However,
these modalities are only approximations of biological reality and come
with their own assumptions and limitations. In time, the inclusion of more
datasets with known ground truth interactions will become available, from
which the limitations and advantages of the different CCC methods will
be better understood.

**This subtask evaluates the methods' ability to predict interactions,
the corresponding of cytokines of which, are inferred to be active in
the target cell types. This subtask focuses
on the prediction of interactions from steady-state, or single-context,
single-cell data.**

## The metrics

Metrics for cell-cell communication aim to characterize how good are
the different scoring methods at prioritizing assumed truth predictions.

* **Odds ratio**: The odds ratio represents the ratio of true and false
positives within a set of prioritized interactions (top ranked hits) versus
the same ratio for the remainder of the interactions. Thus, in this
scenario odds ratios quantify the strength of association between the
ability of methods to prioritize interactions and those interactions
assigned to the positive class.

* **AUPRC**: a single number _[0-1]_ that summarizes the area under the curve where
x is the recall and y is the precision.

## API

### Datasets

Datasets should include cell type annotations in `adata.obs["label"]`, along with some
assumed truth in `adata.uns["ccc_target"]`. The assumed truth could be derived from
various proxies; we refer the reader to [Dimitrov et
al](https://doi.org/10.1038/s41467-022-30755-0) for more details.

`adata.uns["ccc_target"]` should be a Pandas DataFrame containing all the following
columns:

* `response`: `int`, binary response variable _[0; 1]_ indicating whether an interaction is
  assumed to have occurred
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `target`: `str`, name of target cell type in interaction

The datasets should also include a
[NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
in `adata.uns["target_organism"]` - used to convert the (typically human) prior
knowledge of the CCC methods to the corresponding gene homologs.
`adata.X` should contain the raw counts matrix.

### Methods

Methods should predict interactions between cell types without using
`adata.uns["ccc_target"]`. Predicted interactions should be stored in
`adata.uns["ccc_pred"]` as a Pandas DataFrame containing all of the following columns:

* `score`: `float`, score between `-inf` to `+inf` giving a predicted strength of the
  inferred interaction
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `target`: `str`, name of target cell type in interaction

Methods should infer a `score` for each _intersecting interaction_
between a `ligand` and a `target`.
We define _intersecting interactions_ as
those for which the `ligand` genes are present in both the dataset and
the prior-knowledge resource provided by LIANA, while a `target` is any
target cell identity label in the dataset.

The predictions of any method which do not uniquely map
to the columns in `adata.uns["merge_keys"]` are to be **aggregated**.
By default, aggregation is carried as the `max` and `sum`
according to columns in the `merge_keys`.

## Prior-knowledge

Each dataset should be supplemented with a prior knowledge resource of
ligand-receptor interactions, with matching feature IDs.
The resource used in the [Ligand-Target](./cell_cell_communication_ligand_target)
and [Source-Target](./cell_cell_communication_source_target)
tasks was generated as the consensus from multiple manually-curated human
ligand-receptor resources, and includes interactions from
[CellPhoneDB](https://www.nature.com/articles/s41596-020-0292-x),
[CellChatDB](https://www.nature.com/articles/s41467-021-21246-9#disqus_thread),
[ICELLNET](https://www.nature.com/articles/s41467-021-21244-x),
[connectomeDB2020](https://www.nature.com/articles/s41467-020-18873-z),
and [CellTalkDB](https://www.nature.com/articles/s41467-020-18873-z) resources.
All of these were queried via the
[OmniPath database](https://www.embopress.org/doi/full/10.15252/msb.20209923),
and are fixed for each version of LIANA.

To ensure the consistency between the IDs in the dataset and those in the
resource we use a reference map, obtained via BioConductor-v3.15 `org.Hs.eg.db`,
and are provided in `tnbc_wu2021_gene_symbols.csv`.

The prior-knowledge resource is available via the
`cell_cell_communication.utils.ligand_receptor_resource` function, which returns a
DataFrame containing the columns `ligand_genesymbol` and `receptor_genesymbol`, which
correspond to the ligand and receptor genes, respectively. These may contain complexes
with subunits separated with `_`. Hence, **methods should be able to deal with
complex-containing interactions**.

### Metrics

Metrics should evaluate the concordance between `adata.uns["ccc_target"]` and
`adata.uns["ccc_pred"]` to evaluate the success of a method in predicting interactions.
