
# openproblems v2.0.0

A major update to the OpenProblems framework, switching from a Python-based framework to a Viash-based framework. This update includes a complete reorganization of the codebase, with tasks now located in `src/tasks/` and common components in `src/common/`.

Structure:

* src/common: Common components used by all tasks.
* src/datasets: Components for fetching and processing datasets.
* src/tasks/*: Benchmarking tasks
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
* Fix DR baselines (#816)
* set adata.uns['is_baseline'] (#820)
* Copy anndata in metric decorator (#819)
* Don't recompute X_emb and neighborhood graph for baseline datasets (#823)
* Changes in destVI code (#826) (#827)
* Set explicit token permissions (#828)
* Warnings fix (#831)
* Harmonize batch integration dataset APIs (#834)
* new common baselines and cross import (#825)
* jitter baseline patch (#838)
* Add reversed norm order for ALRA in Denoising Task (#835)


# openproblems v0.7.0

## What's Changed
* Fix docker image builds (#758)
* [Dimensionality reduction] Fix normalization in baselines (#760)
* downgrade gtfparse and polars (#766)
* Fix output headers order (#769)
* Convert references to bib (#720)
* fix typo in bibliography path (#774)
* More bibliography typos (#775)
* Pre-normalize dimensionality reduction datasets (#768)
* Add pymde to dimensionality reduction (#767)
* Fix flaky R installations in docker build (#783)
* save initial layer in X for adata_pre (#784)
* Filter datasets by celltype (#770)
* Pass raw counts to neuralee (#779)
* Label projection describe datasets (#776)
* Add missing DR references (#782)
* Bugfix/lowercase GitHub repo owner (#794)
* Upgrade isort (#795)
* Update styler to 1.9.0 (#787)
* [auto] Update docker version (#798)
* Update bslib to 0.4.2 (#759)
* add missing logfc decorator (#796)
* Add ALRA preprocessing identical to literature (#763)
* run CI on PRs only with approving review (#804)
* add new workflow to add status (#805)
* Update bioc/scran to 1.26.2 (#799)
* Specify PR number (#808)
* add magic with reverse norm order (#797)
* Bump pymde from 0.1.15 to 0.1.18 in /docker/openproblems-python-pytorch (#801)
* Update scvi-tools requirement from ~=0.16 to ~=0.19 in /docker/openproblems-r-pytorch (#731)
* Use graph and embedding metrics for feature and embedding subtask (#807)
* Fix typo in dimensionality reduction dataset names (#802)
* add new dataloaders (#792)
* rmse -> distance correlation (#811)
* CPM -> CP10k (#812)
* change multimodal data integration task name to matching modalities  (#778)
* updated scib version (#793)
* Daniel strobl hvg conservation fix (#785)


# openproblems v0.3.1




# openproblems v0.2.1

## Added

* **Pancreas Dataset:** The Pancreas dataset has been added for the "Label Projection" task, including two variants: `pancreas_batch` and `pancreas_random`.
* **MLP Methods:** Multilayer Perceptron (MLP) methods have been introduced for the "Label Projection" task, with variants using scran and log CPM normalization: `mlp_scran` and `mlp_log_cpm`.
* **F1 Score Metric:** The F1 score metric has been added to assess the performance of label projection methods.

### Changed

* **Evaluation Script:** The evaluation script (`evaluate.py`) now sorts results by value in descending order (highest values first) and then by task, dataset, and metric.
* **Normalization in Data Loaders:** Normalization has been removed from the pancreas data loader, and normalization is now handled within the methods themselves.
* **Test Data Size:** The size of the test data for the Zebrafish dataset has been reduced.
* **Dependencies:**  Updated Black formatting and moved some dependencies to optional requirements.
* **Miscellaneous:** Minor code formatting and cleanup.

### Fixed

* **Data Loading:** Addressed a bug in the pancreas data loader to ensure that the test version includes test data.
* **Dataset Checks:**  Enhanced dataset checks to verify the presence of both training and validation data.
* **Normalization:** Ensured that `adata.obsm["mode2"]` is normalized in the MNN methods.

### Removed

* **Dummy and Cheat:** Removed the dummy and cheat metrics and datasets, as they are no longer needed.


# openproblems v0.2

## Added

* **Issue Templates:** Added issue templates to improve bug reporting and streamline the process for proposing new datasets, methods, metrics, and tasks.
* **Zebrafish Dataset:** Added the Zebrafish dataset for the "Label Projection" task, including two variants: `zebrafish_labels` and `zebrafish_random`.

## Changed

* **Dependencies:** Moved `scIB`, `rpy2`, `harmonicalignment`, and `mnnpy` to optional dependencies to reduce installation size for users not using these methods.
* **Travis CI:** The script now explicitly adds the `website/content/results` directory to the git commit.
* **`setup.py`:**  The `extras_require` section has been updated to include the optional dependencies.

## Fixed

* **`TruncatedSVD`:** Fixed issues related to setting the number of components for `TruncatedSVD`.
* **Data Scaling:** Ensured data is scaled appropriately for regression tasks, handling both sparse and dense matrices.
* **Normalization:** Fixed issues related to storing normalization results and normalizing `obsm` data.
* **Harmonic Alignment:** Fixed a bug in the import of `harmonicalignment`.
* **Website:** The site name has been updated in the configuration.
* **Tests:** Improved checks in `test_datasets.py` to ensure datasets have non-zero size.

## Other

* **Version:** The repository version has been bumped to 0.2.0.
* **Miscellaneous:** Minor bug fixes and improvements, including code formatting and handling of macOS-specific files.
* **Documentation:** Documentation related to normalization and issue templates has been updated.

# openproblems v0.1

## Added

* **Website:** A new website built using the Academic Kickstart template has been added to the repository.
* **Results Pages:** Each task now has a dedicated page in the website to display the benchmarking results.
* **Social Sharing:** The website now includes social sharing buttons for easy sharing of content.
* **Data Loader for Zebrafish:** A new data loader for zebrafish data has been added.
* **Normalization Functions:** New normalization functions have been added to `openproblems.tools.normalize`.
* **Logistic Regression Methods:** Logistic regression methods have been added for the "Label Projection" task.
* **Harmonic Alignment Methods:** Harmonic alignment methods have been added for the "Multimodal Data Integration" task.

## Changed

* **Black:** Black configuration now excludes the `website` directory.
* **Tests:** Tests have been updated to reflect changes in data loading and method execution.
* **Version:** The repository version has been bumped to 0.1.
* **`evaluate.py`:**  The script now generates separate Markdown files for each task's results, including frontmatter for the website.

## Fixed

* **Data Loaders:** Issues with the zebrafish data loader and caching have been fixed.
* **`cheat` Method:** The `cheat` method for "Multimodal Data Integration" has been corrected.

## Removed

* **`adata.raw`:** Datasets no longer store raw data in `adata.raw`. Normalization is now handled within methods.


# openproblems v0.0.3

## Added

* **Normalization:**
    * Normalization functions (`log_scran_pooling`, `cpm`, `log_cpm`, `sqrt_cpm`) have been added to `openproblems.tools.normalize`.
    * A `normalizer` decorator has been introduced to streamline the normalization process.
    * Normalization is now performed within methods instead of in the data loaders.

* **Methods:**
    * `harmonic_alignment_sqrt_cpm` and `harmonic_alignment_log_scran_pooling` methods have been added for the "Multimodal Data Integration" task.
    * `logistic_regression_scran` and `logistic_regression_log_cpm` methods have been added for the "Label Projection" task.

* **Testing:**
    * New test modules `test_tools.py` and `test/__init__.py` have been added.
    * Tests now explicitly include their names using `name_func=utils.name_test`.
    * Numba warnings are ignored in tests.
    * The `test_datasets.py` script now checks the size of test datasets.
    * The `test_methods.py` script now verifies that methods operate in-place.

* **Dependencies:**
    * `rpy2`, `scIB`, and `harmonicalignment` have been added to the install requirements.
    * Bioconductor and `scran` are installed in the Travis CI environment.

## Changed

* **Data Handling:**
    * Raw data is no longer checked or normalized in data loaders; normalization is now handled within methods.
    * The `adata.X` is copied in the evaluation script to prevent modifications.
    * The scicar datasets have been consolidated into a single file.
    * The `utils.get_callable_members` function is used to retrieve callable members of modules.

* **Travis CI:**
    * The `evaluate.sh` script has been improved to handle git operations more robustly.

* **Other:**
    * The `cheat` method for "Multimodal Data Integration" now uses PCA with a specified number of components.
    * The repository version has been bumped to 0.0.3.

## Fixed

* A bug in the `harmonic_alignment` import has been fixed.
* A decorator-related issue has been resolved.
* The zebrafish dataset and loader have been removed.
* The `cheat` method for "Label Projection" has been corrected.

## Removed

* The zebrafish dataset and loader have been removed.
* Raw data is no longer stored or checked in datasets.


# openproblems v0.0.2


## Added

* Datasets are now loaded under `openproblems/data` and cached for efficiency.
* A new `dummy` dataset has been added for general use.
* The scicar datasets are now loaded from `openproblems/data` and cached.
* A `loader` decorator has been added to facilitate caching of datasets.

## Changed

* The API has been clarified: methods now load the full dataset even in `test` mode but only return a small subset for faster downstream analysis.
* The scicar datasets have been modified:
    * Genes and peaks with low counts are now removed.
    * The second modality data is stored in `adata.obsm["mode2"]` instead of `adata.uns["mode2"]`.
    * Normalization is done using `scprep` functions.
* The `cheat` method for "Multimodal Data Integration" has been made more memory-efficient.
* The `mnn` and `procrustes` methods for "Multimodal Data Integration" now handle cases where the number of features is smaller than the requested number of SVD components.
* The `mse` metric now converts sparse matrices to dense arrays before computation.
* The `test_datasets.py` script now checks that the test datasets are small enough.
* The `.travis.yml` file has been updated to use `./evaluate.sh` instead of `source evaluate.sh`.
* The minimum required NumPy version has been increased to 1.17.0.
* The `decorator` package has been added to the install requirements.

## Fixed

* A typo in the `evaluate.sh` script has been corrected.
* The path for caching datasets has been fixed.
* A string comparison in the evaluation script has been corrected.
* The `evaluate.sh` script is now explicitly made executable.


# openproblems v0.0.1

First release of OpenProblems.

Tasks:

* Label projection (0 datasets, 0 methods, 1 metric)
* Multimodal data integration (2 datasets, 2 methods, 2 metrics)
