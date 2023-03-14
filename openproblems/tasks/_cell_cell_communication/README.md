# Cell-cell Communication

The growing availability of single-cell data has sparked an increased
interest in the inference of cell-cell communication (CCC),
with an ever-growing number of computational tools developed for this purpose.

Different tools propose distinct preprocessing steps with diverse
scoring functions, that are challenging to compare and evaluate.
Furthermore, each tool typically comes with its own set of prior knowledge.
To harmonize these, [Dimitrov et
al, 2022](https://openproblems.bio/bibliography#dimitrov2022comparison) recently
developed the [LIANA](https://github.com/saezlab/liana) framework, which was used
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

## The subtasks

Subtasks for cell-cell communication are defined by which aspects of a communication
event are detected. Currently, two subtasks aimed at steady-state, or single-context,
CCC predictions are defined:

* [Ligand-Target](./cell_cell_communication_ligand_target): interactions between
  ligand molecules and target cell types; and
* [Source-Target](./cell_cell_communication_source_target): interactions between
  source cell types and target cell types.

More subtasks may be defined that infer communication events on any of the `source`
cell type, the `target` cell type, the `ligand` molecule, and the receptor.
More aspects of the communication may also be added in the future.

## API

### Datasets

Datasets should include cell type annotations in `adata.obs["label"]`, along with some
assumed truth in `adata.uns["ccc_target"]`. The assumed truth could be derived from
various proxies; we refer the reader to [Dimitrov et
al](https://doi.org/10.1038/s41467-022-30755-0) for more details.

`adata.uns["ccc_target"]` should be a Pandas DataFrame containing:

* `response`: `int`, binary response variable _[0; 1]_ indicating whether an interaction
  is assumed to have occurred and at least one of the following columns:

* `source`: `str`, name of source cell type in interaction
* `target`: `str`, name of target cell type in interaction
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `receptor`: `str`, gene symbol of the receptor in an interaction

The datasets should also include a [NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
in `adata.uns["target_organism"]` - used to convert the (typically human) prior
knowledge of the CCC methods to the corresponding gene homologs.
`adata.X` should contain the raw counts matrix.

### Methods

Methods should predict interactions between cell types without using
`adata.uns["ccc_target"]`. Predicted interactions should be stored in
`adata.uns["ccc_pred"]` as a Pandas DataFrame containing:

* `score`: `float`, score between `-inf` to `+inf` giving a predicted strength of the
  inferred interaction

and at least two of the following columns:

* `source`: `str`, name of source cell type in interaction
* `target`: `str`, name of target cell type in interaction
* `ligand`: `str`, gene symbol of the ligand in an interaction
* `receptor`: `str`, gene symbol of the receptor in an interaction

The relevance of these columns is determined by the subtask in question
via `adata.uns["merge_keys"]`, a list of at least two columns from the
aforementioned columns corresponding to the assumed
truth in `adata.uns["ccc_target"]`.

Methods should infer a score for each _intersecting interaction_,
where these represent the intersecting columns between `adata.uns["ccc_pred"]` and
`adata.uns["ccc_target"]`.

In case, `ligand` and/or `receptor` columns are present
in `adata.uns["ccc_target"]`, we further define _intersecting interactions_ as
those for which the relevant genes are present in both the dataset and
the prior-knowledge resource provided by LIANA.

The predictions of any method which do not uniquely map
to the columns in `adata.uns["merge_keys"]` are to be **aggregated**.
By default, aggregation is carried as the `max` and `sum`
according to columns in the `merge_keys`.

The prior-knowledge resource is available via the
`cell_cell_communication.utils.ligand_receptor_resource` function, which returns a
DataFrame containing the columns `ligand_genesymbol` and `receptor_genesymbol`, which
correspond to the ligand and receptor genes, respectively. These may contain complexes
with subunits separated with `_`. Hence, **methods should be able to deal with
complex-containing interactions**.

### Prior-knowledge Resource

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

### Metrics

Metrics should evaluate the concordance between `adata.uns["ccc_target"]` and
`adata.uns["ccc_pred"]` to evaluate the success of a method in predicting interactions.

### Examples

The triple [negative breast cancer dataset](
https://www.nature.com/articles/s41588-021-00911-1) (`tnbc_wu2021`) portrays
benchmark truth in the form of inferred cytokine activities in the target cell
types, as such in addition to the `response` column `adata.uns["ccc_target"]`,
also contains `ligand` and `target` columns.

In the case of the [murine brain dataset](
https://www.nature.com/articles/nn.4216) (`allen_brain_atlas`), we assume that
spatially-adjacent cell types are more likely to interact, hence interactions
between them should be preferentially detected.
Consequently, `adata.uns["ccc_target"]` contains `source`, `target`,
and `response` columns, but no `ligand` or `receptor` columns.
