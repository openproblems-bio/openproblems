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
[LIANA](https://github.com/saezlab/liana) framework, which is used
as a foundation for this task.

The challenges in evaluating the tools are further exacerbated by the
lack of an appropriate gold standard to benchmark the performance of
CCC methods. To solve this, Dimitrov et al use alternative data modalities,
including the spatial proximity of cell types and inferred downstream
cytokine activities, to generate an inferred ground truth. However,
these modalities are only approximations of biological reality and come
with their own assumptions and limitations. In time, the inclusion of more
datasets with known ground truth interactions will become available, from
which the limitations and advantages of the different CCC methods will
be better understood.

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

Datasets should include cell type annotations in `adata.obs["label"]`,
along with some assumed truth in `adata.uns["ccc_target"]` with a binary
`response` vector [0; 1], and the corresponding columns relevant for the
benchmark. The assumed truth could essentially be in
different forms, and we refer the reader to [our recent publication for more
details](https://rdcu.be/cSs92).

For example, the triple negative breast cancer dataset (`tnbc_wu2021`) portrays
benchmark truth in the form of inferred cytokine activities in the target cell
types, as such in addition to the `response` column `adata.uns["ccc_target"]`,
also contains `ligand` and `target` columns with which we can join the assumed
truth to the output CCC predictions in `adata.uns['ccc']`.
In the case of the murine brain dataset (`allen_brain_atlas`), we assume that
spatially-adjacent cell types are more likely to interact, hence interactions
between them should be preferentially detected.
Consequently, `adata.uns["ccc_target"]` contains `source`, `target`,
and `response`columns.

The datasets should also include a
[NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
in `adata.uns["target_organism"]` - used to convert the (typically human) prior
knowledge of the CCC methods to the corresponding gene homologs.
`adata.X` should contain the raw counts matrix.

`adata.uns["ccc_target"]` should be a Pandas DataFrame containing the following
columns:
- `response`: int, binary response variable indicating whether an interaction
is assumed to have occurred
- `source`: str, name of source cell type in interaction
- `target`: str, name of target cell type in interaction
- `ligand`: str, gene symbol of the ligand in an interaction
- `receptor`: str, gene symbol of the receptor in an interaction (may be null)
- `score`: float, score between -inf to +inf giving an inferred strength of
the inferred interaction

Metrics should evaluate the ability of methods to preferentially
detect assumed truth CCC events.
