
# openproblems v2.0.0

A major update to the OpenProblems framework, switching from a Python-based framework to a Viash-based framework. This update includes a complete reorganization of the codebase, with tasks now located in `src/tasks/` and common components in `src/common/`.

Structure:

* src/common: Common components used by all tasks.
* src/datasets: Components for fetching and processing datasets.
* src/tasks/*: Benchmarking tasks
* src/migration: Components for migrating from v1 to v2.
* src/wf_utils: Workflow utilities.

Migrates tasks:

* batch_integration
* denoising
* dimensionality_reduction
* match_modalities
* predict_modality
* spatial_decomposition
* spatially_variable_genes

# openproblems v0.8.0

## What's Changed
* Fix DR baselines by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/816
* set adata.uns['is_baseline'] by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/820
* Copy anndata in metric decorator by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/819
* Don't recompute X_emb and neighborhood graph for baseline datasets by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/823
* Changes in destVI code (#826) by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/827
* Set explicit token permissions by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/828
* Warnings fix by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/831
* Harmonize batch integration dataset APIs by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/834
* new common baselines and cross import by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/825
* jitter baseline patch by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/838
* Add reversed norm order for ALRA in Denoising Task by @wes-lewis in https://github.com/openproblems-bio/openproblems/pull/835


# openproblems v0.7.0

## What's Changed
* Fix docker image builds by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/758
* [Dimensionality reduction] Fix normalization in baselines by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/760
* downgrade gtfparse and polars by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/766
* Fix output headers order by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/769
* Convert references to bib by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/720
* fix typo in bibliography path by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/774
* More bibliography typos by @scottgigante in https://github.com/openproblems-bio/openproblems/pull/775
* Pre-normalize dimensionality reduction datasets by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/768
* Add pymde to dimensionality reduction by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/767
* Fix flaky R installations in docker build by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/783
* save initial layer in X for adata_pre by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/784
* Filter datasets by celltype by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/770
* Pass raw counts to neuralee by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/779
* Label projection describe datasets by @mxposed in https://github.com/openproblems-bio/openproblems/pull/776
* Add missing DR references by @rcannood in https://github.com/openproblems-bio/openproblems/pull/782
* Bugfix/lowercase GitHub repo owner by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/794
* Upgrade isort by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/795
* Update styler to 1.9.0 by @github-actions in https://github.com/openproblems-bio/openproblems/pull/787
* [auto] Update docker version by @github-actions in https://github.com/openproblems-bio/openproblems/pull/798
* Update bslib to 0.4.2 by @github-actions in https://github.com/openproblems-bio/openproblems/pull/759
* add missing logfc decorator by @dbdimitrov in https://github.com/openproblems-bio/openproblems/pull/796
* Add ALRA preprocessing identical to literature by @wes-lewis in https://github.com/openproblems-bio/openproblems/pull/763
* run CI on PRs only with approving review by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/804
* add new workflow to add status by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/805
* Update bioc/scran to 1.26.2 by @github-actions in https://github.com/openproblems-bio/openproblems/pull/799
* Specify PR number by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/808
* add magic with reverse norm order by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/797
* Bump pymde from 0.1.15 to 0.1.18 in /docker/openproblems-python-pytorch by @dependabot in https://github.com/openproblems-bio/openproblems/pull/801
* Update scvi-tools requirement from ~=0.16 to ~=0.19 in /docker/openproblems-r-pytorch by @dependabot in https://github.com/openproblems-bio/openproblems/pull/731
* Use graph and embedding metrics for feature and embedding subtask by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/807
* Fix typo in dimensionality reduction dataset names by @lazappi in https://github.com/openproblems-bio/openproblems/pull/802
* add new dataloaders by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/792
* rmse -> distance correlation by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/811
* CPM -> CP10k by @scottgigante-immunai in https://github.com/openproblems-bio/openproblems/pull/812
* change multimodal data integration task name to matching modalities  by @LuckyMD in https://github.com/openproblems-bio/openproblems/pull/778
* updated scib version by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/793
* Daniel strobl hvg conservation fix by @danielStrobl in https://github.com/openproblems-bio/openproblems/pull/785
