# openproblems v2.1.0

## Minor changes

- Add the CELLxGENE immune cell atlas dataset as a common test resource (PR #907)

# openproblems v2.0.0

A major update to the OpenProblems framework, switching from a Python-based framework to a Viash + Nextflow-based framework. This update features the same concepts as the previous version, but with a new implementation that is more flexible, scalable, and maintainable.

Most relevant parts of the overall structure:

* `src/tasks`: Benchmarking tasks:
  - `batch_integration`: Batch integration
  - `denoising`: Denoising
  - `dimensionality_reduction`: Dimensionality reduction
  - `match_modalities`: Match modalities
  - `predict_modality`: Predict modality
  - `spatial_decomposition`: Spatial decomposition
  - `spatially_variable_genes`: Spatially variable genes

* `src/datasets`: Components for creating common datasets. Loaders:
  - `cellxgene_census`: Query cells from a CellxGene Census
  - `openproblems_neurips2021_bmmc`: Fetch a dataset from the OpenProblems NeurIPS2021 competition
  - `openproblems_neurips2022_pbmc`: Fetch a dataset from the OpenProblems NeurIPS2022 competition
  - `openproblems_v1`: Fetch a legacy OpenProblems v1 dataset
  - `openproblems_v1_multimodal`: Fetch a legacy OpenProblems v1 multimodal dataset
  - `tenx_vision`: Fetch a and convert 10x Visium dataset
  - `zenodo_spatial`: Fetch and process an Anndata file containing DBiT seq, MERFISH, seqFISH, Slide-seq v2, STARmap, and Stereo-seq data from Zenodo.
  - `zenodo_spatial_slidetags`: Download a compressed file containing gene expression matrix and spatial locations from zenodo.

* `src/common`: Common components used by all tasks.
  - `check_dataset_schema`: Check whether an h5ad dataset adheres to a dataset schema
  - `check_yaml_schema`: Check whether a YAML adheres to a JSON schema
  - `comp_tests`: Reusable component unit tests
  - `create_component`: Create a component Viash component.
  - `create_task_readme`: Create a README for an OpenProblems task.
  - `extract_metadata`: Extract the `.uns` metadata from an h5ad file.
  - `helper_functions`: Commonly used helper functions in Python or in R,
  - `process_task_results`: Process the raw tasks results (containing raw logs, unprocessed component configs, and various metrics) into nicely formatted task results.
  - `schemas`: JSON schemas for YAML files in the repository
  - `sync_test_resources`: Synchronise the test resources from s3 to resources_test

For more information related to the structure of this repository, see the [documentation](https://openproblems.bio/documentation/reference/openproblems/).


# openproblems v1.0.0

Note: This changelog was automatically generated from the git log.

## New functionality
- Added `cell2location` to the `spatial_decomposition` task.
- Added nearest-neighbor ranking matrix computation to `_utils`.
- Datasets now store nearest-neighbor ranking matrix in `adata.obsm["X_ranking"]`.
- Added support for parsing Nextflow output and generating benchmark results for the website.
- Added `max_samples` parameter to `qlocal`, `qglobal`, `qnn_auc`, `lcmc`, `qnn`, and `continuity` metrics to allow for subsampling of data for faster computation.
- Added new scArches based methods: `scarches_scanvi_xgb_all_genes` and `scarches_scanvi_xgb_hvg`.
- Added `prediction_method` parameter to `_scanvi_scarches` to specify prediction method.
- Added `_pred_xgb` function to perform XGBoost prediction based on latent representations.
- Added `obsm` parameter to `_xgboost` function to allow specifying the embedding space for XGBoost training.

## Major changes
- Updated `scvi-tools` to version `0.20` in both Python and R environments.
- Updated datasets to include nearest-neighbor ranking matrix.
- Modified dimensionality reduction task to include nearest-neighbor ranking matrix computation in dataset generation.
- The website update workflow was refactored to use a new workflow using json instead of markdown.
- Updated the website generation process to remove duplicate BibTex entries.
- Added a new `parse_metadata.py` script for generating metadata for the website.
- Added a new function to `openproblems.utils.py` to get the member ID of a task, dataset, method or metric.
- Removed the redundant computation and storage of the nearest-neighbor ranking matrix in datasets.

## Minor changes
- Updated method names to be shorter and more consistent across tasks.
- Improved method summaries for clarity.
- Updated JAX and JAXlib versions to 0.4.6.
- Updated dependencies to support new versions of Snakemake and GitPython.
- Removed code related to "nbt2022-reproducibility" repo and merged it into the main website.
- Updated the schema for benchmark results to include submission time, code version, and resource usage metrics.
- Improved error handling and added logging to the parsing script.
- Removed the "raw.json" file from the results directory and merged all data into a single "results.json" file.
- Updated the workflow to upload the final results to the website's results directory instead of the data directory.
- Removed unnecessary code and refactored the parsing script for better readability.
- Added unit tests for the new parsing script.
- Updated the `run_tests` workflow to skip testing on the `test_website` branch.
- Updated the `run_tests` workflow to skip testing on the `test_process` branch.
- Updated the `create-pull-request` step to set the author for the pull request.
- Updated the `run_tests` workflow to skip testing on pull request reviews.
- Updated the `update_website_content` workflow to update the website on the `main` branch.
- Updated the `main.bib` file to fix a typo.
- Removed extraneous headings from task README files.
- Updated `generate_test_matrix.py` to use the new `openproblems.utils.get_member_id` function.
- Updated the website generation process to copy BibTex files to the correct location.
- Updated the `process_requires` section in `setup.py` to include `gitpython`.
- Updated git commit hash generation for openproblems functions.
- Modified `_xgboost` to allow for specifying `tree_method`.
- Modified `_scanvi_scarches` to consistently use `unlabeled_category`.
- Modified `_scanvi_scarches` to remove unnecessary copying of `labels`.
- Removed `_scanvi_scarches` functions that were redundant with `_scanvi_scarches`.
- Removed unused `_scanvi` functions.
- Modified `_scanvi_scarches` to allow for specifying `prediction_method` and handle `unlabeled_category` consistently.

## Documentation
- Improved the documentation of the `auprc` metric.
- Improved the documentation of the `cell2location` methods.
- Document sub-stub task behaviour

## Bug fixes
- Fixed an error in `neuralee_default` where the `subsample_genes` argument could be too small.
- Fixed an error in `knn_naive` where the `is_baseline` argument was set to `False`.
- Fixed calculation of ranking matrix in `_utils` to include ties.
- Fixed a bug in `load_tenx_5k_pbmc()` where a warning about non-unique variable names was being raised.
- Removed the unused `_utils.py` file.
- Removed the `X_ranking` entry from the `obsm` attribute of datasets.
- The `_fit()` function in `nn_ranking.py` now subsamples the data if `max_samples` is specified.
- The `nn_ranking` metrics now use subsampling in the `_fit()` function to improve performance.
- Fixed the git hash generation for openproblems functions
- Fixed a warning about `pkg_resources` being deprecated
- Removed unnecessary `fetch-depth: 1` from workflow
- Fixed potential issue in `_scanvi_scarches` where `labels_pred` could be overwritten
- Fixed potential issue in `_pred_xgb` where `num_round` wasn't being used correctly
- Fixed an issue where baseline methods were not being filtered correctly from the benchmark results.
- Fixed an issue where metrics with all NaN values were not being removed from the benchmark results.
- Fixed an issue where some metrics were not being parsed correctly from the Nextflow output.
- Fixed an issue where the "mean_score" field was not being calculated correctly for each method.
- Fixed an issue where the "code_version" field was not being populated correctly for each method.
- Fixed an issue where the "submission_time" field was not being populated correctly for each method.
- Fixed an issue where the resource usage metrics were not being parsed correctly from the Nextflow output.
- Updated the `run_tests` workflow to skip testing on the `test_website` branch.
- Updated the `run_tests` workflow to skip testing on the `test_process` branch.
- Updated the `create-pull-request` step to set the author for the pull request.
- Updated the `run_tests` workflow to skip testing on pull request reviews.
- Updated the `update_website_


# openproblems v0.8.0

Note: This changelog was automatically generated from the git log.

## New functionality
- Added the zebrafish_labs dataset to the dimensionality reduction task.
- Added the `diffusion_map` method to the dimensionality reduction task.
- Added the `spectral_features` method to the dimensionality reduction task, which uses diffusion maps to create embedding features.
- Added the `distance_correlation_spectral` metric to the dimensionality reduction task, which evaluates the similarity of the high-dimensional Laplacian eigenmaps on the full data matrix and the dimensionally-reduced matrix.
- Added baseline methods for batch integration: no integration, random integration, random integration by cell type, random integration by batch.
- Added `alra_sqrt_reversenorm`, `alra_log_reversenorm` methods for ALRA with reversed normalization order.
- Added `celltype_random_embedding_jitter` method to randomize embedding with jitter.

## Minor changes
- Improved the `density_preservation` metric calculation.
- Updated the `distance_correlation` metric to use the new `diffusion_map` method.
- Increased the default number of components used for `distance_correlation_spectral` to 1000.
- Made metrics more robust by copying the AnnData object before passing it to the metric function.
- Added `is_baseline` flag to `adata.uns` in `method` decorator.
- Added `is_baseline` field to `adata.uns` for all methods.
- Increased default values for `max_epochs_sp` and `max_epochs_sc` in `destvi` method.
- Changed default value of `early_stopping_monitor` to `elbo_validation` from `reconstruction_loss_train` in `destvi` method.
- Added `train_size` and `validation_size` arguments to the `sc_model.train` call in `destvi` method.
- Added `batch_size` and `plan_kwargs` arguments to the `st_model.train` call in `destvi` method.
- Refactor ALRA methods for improved clarity and consistency.
- Added tests for new ALRA methods with reversed normalization order.
- Added jitter parameter to `_random_embedding` function.
- Updated `celltype_random_embedding` to use `jitter=None` in `_random_embedding`.
- Removed unnecessary parameters from the `sample_dataset` function.
- Removed unnecessary checks for PCA and neighbors in the `check_dataset` function.
- Updated `pytest.ini` to ignore deprecation warning related to `pkg_resources`.
- Added permission to all workflows to read and write contents
- Added permission to write pull requests to several workflows
- Added permission to write packages to the `run_tests` workflow.

## Bug fixes
- Fixed a bug in `density_preservation` that caused it to return 0 when there were NaN values in the embedding.
- Removed unused `true_features_log_cp10k` and `true_features_log_cp10k_hvg` methods.
- Removed unnecessary imports in metrics.
- Removed unnecessary `neighbors` calls in metrics.
- Removed unused `_get_split` function.
- Added `embedding_to_graph` and `feature_to_graph` functions for graph-based metrics.
- Added `get_split` function for metrics that require splitting data into training and testing sets.
- Added `feature_to_embedding` function for embedding-based metrics.
- Fixed issue where baseline methods were not properly documented.
* Increased default maximum epochs for spatial models to improve performance.
* Improved training parameters for both spatial and single-cell models to improve stability and performance.
* Updated validation metric used for early stopping in spatial model to improve training quality.

## Documentation
- Updated documentation to clarify that the AnnData object passed to metric functions is a copy.
- Updated the documentation for batch integration tasks to reflect the change in the expected format of the dataset objects.

## Major changes
- Moved baseline methods from individual task modules to a common module.
- Removed redundant baseline methods from individual task modules.
- Increased default values for `max_epochs_sp` and `max_epochs_sc` in `destvi` method.


# openproblems v0.7.4

Note: This changelog was automatically generated from the git log.

## New functionality
- Added metadata for all datasets, methods, and metrics.

## Major changes
- Updated nf-openproblems to v1.10.

## Minor changes
- Added a new `docker_pull` rule to the Snakemake workflow to pull Docker images.
- Added a new `docker` rule to the Snakemake workflow to build Docker images.
- Changed the `pytest` command to include coverage for the `test` directory.
- Added new environment variables for the TOWER_TEST_ACTION_ID and TOWER_FULL_ACTION_ID to the Snakemake workflow.
- Updated the `scripts/install_renv.R` script to increase the number of retry attempts.


# openproblems v0.7.3

Note: This changelog was automatically generated from the git log.

## Minor changes
- Updated `scib` version to `1.1.3` in `docker/openproblems-r-extras/requirements.txt` and `docker/openproblems-r-pytorch/requirements.txt`.
## Bug fixes
- Added `pytest-timestamper` to test dependencies for better debugging.


# openproblems v0.7.2

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed an issue where pymde did not work on sparse data.


# openproblems v0.7.1

Note: This changelog was automatically generated from the git log.

## Minor changes
- Added `hvg_unint` and `n_genes_pre` to the lung batch.


# openproblems v0.7.0

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a bibtex file `main.bib` for storing all references cited in the repository.
- Added a section on adding paper references to `CONTRIBUTING.md` explaining how to add entries to `main.bib` and link to them in markdown documents.
- Added new baseline methods for dimensionality reduction: "True Features (logCPM)", "True Features (logCPM, 1kHVG)".
- Added `alra_log` method, which implements ALRA with log normalization.
- Added `alra_sqrt` method, which implements ALRA with square root normalization.
- Added PyMDE dimensionality reduction methods
- Added citations for Chen et al. (2009) "Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis", Kraemer et al. (2018) "dimRed and coRanking - Unifying Dimensionality Reduction in R", Lee et al. (2009) "Quality assessment of dimensionality reduction: Rank-based criteria", Lueks et al. (2011) "How to Evaluate Dimensionality Reduction? - Improving the Co-ranking Matrix", Szubert et al. (2019) "Structure-preserving visualisation of high dimensional single-cell datasets", and Venna et al. (2006) "Local multidimensional scaling".
- Added `install_renv.R` script to install R packages using `renv` with retries
- Added a new metric to evaluate the conservation of highly variable genes (HVGs) after batch integration.
- Added support for lung data from Vieira Braga et al.
- Added `magic_reverse_norm` and `magic_approx_reverse_norm` methods which reverse the order of normalization and transformation in the MAGIC algorithm.
- Added a new workflow to comment on pull request status.

## Major changes
- Updated the `openproblems` repository to cite papers using bibtex references.
- Renamed `alra` method to `alra_sqrt`.
- Updated `spacexr` to latest version.
- Added `fc_cutoff` and `fc_cutoff_reg` parameters to `rctd` method to control minimum log-fold-change for genes in the normalization and RCTD steps.
- Renamed the "multimodal_data_integration" task to "matching_modalities".
- Bumped version to 0.7.0.

## Minor changes
- Added BibTex references to all data loaders in `openproblems/data`.
- Added BibTex references to all methods in `openproblems/tasks`.
- Added BibTex references to all metrics in `openproblems/tasks`.
- Updated `update_website_content.yml` to copy `main.bib` to the Open Problems website.
- Added a BibTeX Tidy hook to `.pre-commit-config.yaml`.
- Updated `scvi-tools` version to `~0.19` in both `openproblems-python-pytorch` and `openproblems-r-pytorch` dockerfiles.
- Updated `cell2location` version to `47c8d6dc90dd3f1ab639861e8617c6ef0b62bb89` in the `openproblems-python-pytorch` dockerfile.
- Updated `bslib` to version 0.4.2.
- Updated `htmltools` to version 0.5.4.
- Updated the `alra_sqrt` method to use square root normalization.
- Updated the `alra_log` method to use log normalization.
- Updated the method names to reflect the normalization used.
- Updated dependencies for `gtfparse` and `polars`.
- Added PyMDE dependency to requirements.txt
- Updated the API to specify that datasets should provide log CPM-normalized counts in `adata.X


# openproblems v0.6.1

Note: This changelog was automatically generated from the git log.

## New functionality
- Added `cell2location_detection_alpha_1` method, which uses `detection_alpha=1` and a hard-coded reference.
- Added a new parameter `hard_coded_reference` to `cell2location_detection_alpha_1` method.
- Added a new baseline method for dimensionality reduction using high-dimensional Laplacian Eigenmaps.
- Added organism metadata to datasets.
- Added a new image, `openproblems-python-bedtools`, to contain packages required for running `pybedtools` and `pyensembl` Python packages.
- Added support for TensorFlow 2.9.0.
- Added a new schema for storing results in JSON format.
- Added a new function to parse Nextflow trace files to this JSON schema.
- Added `rmse_spectral` metric, which calculates the root mean squared error (RMSE) between high-dimensional Laplacian eigenmaps on the full (or processed) data matrix and the dimensionally-reduced matrix.
- Added new methods to LIANA: `magnitude_max`, `magnitude_sum`, `specificity_max`, and `specificity_sum`.
- Added `aggregate_how` parameter to `liana` R function to allow aggregation by "magnitude" or "specificity".
- Added `top_prop` parameter to `odds_ratio` metric to allow specifying the proportion of interactions to consider for calculating the odds ratio.

## Major changes
- Removed unused `openproblems-python-batch-integration` docker image.
- Moved `scanorama`, `bbknn`, `scVI`, `mnnpy` and `scib` from `openproblems-python-batch-integration` to `openproblems-r-pytorch`.
- Moved `cell2location`, `molecular-cross-validation`, `neuralee`, `tangram` and `phate` from `openproblems-python-extras` to `openproblems-python-pytorch`.
- Moved `pybedtools`, `pyensembl` and `scalex` from `openproblems-python-extras` to `openproblems-python-pytorch`.
- Moved `dca` and `keras` from `openproblems-python-tf2.4` to `openproblems-python-tensorflow`.
- Added `openproblems-python-bedtools` docker image.
- Added `openproblems-python-tensorflow` docker image.
- Added `openproblems-python-pytorch` docker image.
- Moved `harmony-pytorch` from `openproblems-r-extras` to `openproblems-r-pytorch`.
- Added `openproblems-r-pytorch` docker image.
- Updated `anndata2ri` version in `openproblems-r-base`.
- Updated `kBET` version in `openproblems-r-extras`.
- Updated `scib` version in `openproblems-r-extras`.
- Updated `scvi-tools` version in `openproblems-r-pytorch`.
- Updated `torch` version in `openproblems-r-pytorch`.
- Moved the `codecov` action to run only on success
- Updated the workflow to upload coverage reports to GitHub Actions as an artifact
- Renamed the `run_benchmark` job to `setup_benchmark`.
- Added a new `run_benchmark` job that runs after `setup_benchmark`.
- Moved the benchmark running logic from the `run_benchmark` job to the new `run_benchmark` job.
- Added a `setup-environment` step to `setup_benchmark` job.
- Added outputs to the `setup_benchmark` job.
- Renamed the `nbt2022-reproducibility` to `website-experimental`

## Minor changes
- Updated `numpy` and `scipy` dependencies in setup.py.
- Updated `scikit-learn`, `louvain`, `python-igraph`, `decorator` and `colorama` dependencies in setup.py.
- Improved Docker image caching.
- Removed the `counts` layer from the `immune_cells`, `pancreas` datasets, and the `batch_integration_feature` task.
- Removed the `counts` layer from `generate_synthetic_dataset` functions in spatial decomposition datasets.
- Updated the `normalize` functions to not modify the data in place.
- Updated the `log_cpm_hvg` function to annotate HVGs instead of subsetting the data.
- Updated the `_high_dim` function in the `nn_ranking` metric to subset to HVGs.
- Updated the `dimensionality_reduction` task README to clarify the role of the `highly_variable` key.
- Reduced the random noise added to the one-hot embedding in the `_random_embedding` function from (-0.1, 0.1) to (-0.01, 0.01).
- Removed `high_dim_pca` and `high_dim_spectral` methods.
- Updated the `random_features` method to use the `check_version` function.
- Moved raw output files from website to the NBT 2022 reproducibility repository.
- Updated the `process_results.yml` workflow to include the NBT 2022 reproducibility repository.
- Updated the `run_tests.yml` workflow to skip tests when pushing to specific branches.
- Removed `# ci skip` from commit message in CI workflow.
- Removed redundant file deletion from `process_results.yml` workflow.
- Added `update_website_content.yml` workflow to update benchmark content on the website repository.
- Modified the `process_results.yml` workflow to update website content based on results.
- Changed the `update_website_content.yml` workflow to trigger on both the `main` and `test_website` branches.
- Updated workflow to push changes to the website only if there are changes to the website content.
- Added environment variable to track changes.
- Removed unused git command.
- Decreased number of samples for testing.
- Updated `igraph` to 0.10.* in `setup.py`.
- Updated `anndata2ri` to 1.1.* in `openproblems-r-base/README.md`.
- Updated `kBET` to `a10ffea` in `openproblems-r-extras/r_requirements.txt`.
- Updated `scib` to `f0be826` in `openproblems-r-extras/requirements.txt`.
- Updated `harmony-pytorch` to 0.1.* in `openproblems-r-pytorch/requirements.txt`.
- Updated `torch` to 1.13.* in `openproblems-r-pytorch/requirements.txt`.
- Updated `scanorama` to 1.7.0 in `openproblems-r-pytorch/requirements.txt`.
- Updated `scvi-tools` to 0.16.* in `openproblems-r-pytorch/requirements.txt`.
- Updated the `regulatory_effect_prediction` task to use


# openproblems v0.6.0

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a new dataset: "Pancreas (inDrop)"
- Added a new function: "pancreas"
- Added a new utility function: "utils.split_data"
- Added `tabula_muris_senis_lung_random` dataset.
- Added `celltype_random_embedding` baseline method for batch integration embedding.
- Added `celltype_random_graph` baseline method for batch integration graph.
- Added a new argument `sctransform_n_cells` to the seuratv3 function to allow users to specify the number of cells used to build the negative binomial regression in the SCTransform function.
- Added a new sample dataset that is smaller and more efficient than the previous one.
- Added a "mean score" metric to the results table.
- Added support for loading the sample dataset in `load_sample_data`.
- Added support for running benchmarks on pull requests.
- Added a new workflow for creating a test matrix.
- Added a new script to generate a test matrix for the `run_tester` workflow.
- Added a new script for cleaning up runner diskspace.
- Added support for uploading docker images to ECR.

## Minor changes
- Added `tabula_muris_senis` dataset to `openproblems/tasks/denoising/datasets/__init__.py`.
- Updated `styler` to version 1.8.1.
- Updated the method for normalizing scores to correctly account for baseline method scores.
- Improved the way NaN and infinite values are handled in the ranking calculation.
- Removed redundant code that was previously used to upload results and markdown artifacts to test.
- Removed the raw output files from the website data directory.
- Updated the list of reviewers for the pull request to include more relevant team members.
- Changed the reference to "Code" to "Library" in the JSON output to better reflect the data presented.
- Added a check to ensure that the task has a minimum number of non-baseline methods before processing results.
- Removed the check to ensure that the task has a minimum number of methods before processing results.
- Removed redundant code that was previously used to handle incomplete tasks.
- Updated the workflow to use a consistent version of Python across all jobs.
- Updated flake8 dependency to `https://github.com/pycqa/flake8`.
- Improved random embedding for `celltype_random_embedding` and `celltype_random_graph`.
- Removed `pip check` from Dockerfile.
- Updated code to use a more consistent random number generator.
- Updated liana code to inverse the distribution of the aggregate rank.
- Improved the logic in `odds_ratio` to ensure that the numerator/denominator is not zero.
- Removed unnecessary NXF_DEFAULT_DSL from `run_tester` workflow.
- Increased the number of cells used to build the negative binomial regression in the SCTransform function from 3000 to 5000.
- Adjusted the default values for `n_pca` and `sctransform_n_cells` in the seuratv3 function for test and non-test cases.
- Updated the seuratv3_wrapper.R script to pass the `sctransform_n_cells` argument to the SCTransform function.
- Moved the sample dataset from the `multimodal` folder to the `sample` folder.
- Refactored the sample data generation to be more efficient.
- Modified the `compute_ranking` function to calculate and add the "mean score" to the `dataset_results` dictionary.
- Updated the `dataset_results_to_json` function to include the "mean score" in the results table.
- Updated the pull request template to reflect recent changes and improvements in the workflow.
- Updated the workflow to include a new `test_full_benchmark` branch.
- Removed redundant code from the workflow.


# openproblems v0.5.21

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a new metric, AUPRC, for evaluating cell-cell communication predictions.
- Added support for aggregating method scores using "max" and "sum" operations.
- Implemented a new method, true events, which predicts all possible interactions.
- Added a new method, random events, which randomly predicts interactions.
- Implemented LIANA, CellPhoneDB, Connectome, Log2FC, NATMI, and SingleCellSignalR methods with the option to aggregate scores using "max" or "sum."
- Added LIANA, CellPhoneDB, Connectome, Log2FC, NATMI, and SingleCellSignalR methods to the cell-cell communication ligand-target task.
- Added LIANA, CellPhoneDB, Connectome, Log2FC, NATMI, and SingleCellSignalR methods to the cell-cell communication source-target task.

## Bug fixes
- Fixed a bug where the odds ratio metric was not handling cases where the numerator or denominator was zero.

## Minor changes
- Updated the IRkernel package version in the R base docker image to 1.3.1.
- Updated the saezlab/liana package version in the R extras docker image to 0.1.7.
- Updated the boto3 package version in the main docker image to 1.26.*.
- Added a check to the cell-cell communication dataset validation to ensure that there are no duplicate entries in the target data.
- Updated the documentation for the cell-cell communication ligand-target task.
- Updated the documentation for the cell-cell communication source-target task.


# openproblems v0.5.20

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed an issue where a sparse matrix was not being converted to CSR format.
- Fixed a bug in `docker_run.sh` where pip check was not being executed.

## Minor changes
- Updated `pkgload` to version 1.3.1.


# openproblems v0.5.19

Note: This changelog was automatically generated from the git log.

## Minor changes
- Converted sparse matrix to csr format.


# openproblems v0.5.18

Note: This changelog was automatically generated from the git log.

## Minor changes
- Converted sparse matrices to CSR format.


# openproblems v0.5.17

Note: This changelog was automatically generated from the git log.




# openproblems v0.5.16

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed a bug where the bioconductor version was incorrect.
- Fixed a bug where the matrix in obs was incorrect.
## Minor changes
- Updated the scran package to version 1.24.1.
- Updated the batchelor and scuttle packages.


# openproblems v0.5.15

Note: This changelog was automatically generated from the git log.




# openproblems v0.5.14

Note: This changelog was automatically generated from the git log.

## Major changes
- Updated workflow to run tests against `prod` branch.


# openproblems v0.5.13

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Skip benchmark if tester fails.


# openproblems v0.5.12

Note: This changelog was automatically generated from the git log.

## New functionality
- Explicitly push prod images on tag

## Documentation
- Added short metric descriptions to README

## Minor changes
- Added labels tests


# openproblems v0.5.11

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Reverted bump of louvain to 0.8, which caused issues.

## Minor changes
- Updated torch requirement to 1.13 in the openproblems-r-pytorch docker.


# openproblems v0.5.10

Note: This changelog was automatically generated from the git log.

## New functionality
- Added support for SCALEX version 1.0.2.

## Minor changes
- Updated RcppAnnoy to version 0.0.20.
- Updated SageMaker requirement to version 2.116.*.

## Bug fixes
- Fixed a bug in the `docker_hash` function, which now returns a string instead of an integer.
- Fixed a bug in the `scalex` method, which now correctly handles the `outdir` parameter.


# openproblems v0.5.9

Note: This changelog was automatically generated from the git log.

## Minor changes
- Update rpy2 requirement from <3.5.5 to <3.5.6
- Update ragg to 1.2.4
## Bug fixes
- Don't fail job if hash fails


# openproblems v0.5.8

Note: This changelog was automatically generated from the git log.

## Minor changes
- Updated scIB to 77ab015.


# openproblems v0.5.7

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a new batch integration subtask for corrected feature matrices.
- Added a new sub-task for batch integration, "batch integration embed", which includes all methods that output a joint embedding of cells across batches.
- Added a new sub-task for batch integration, "batch integration graph", which includes all methods that output a cell-cell similarity graph (e.g., a kNN graph).

# openproblems v0.5.6

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed an issue where the `::` in branch names would cause problems.
- Fixed an issue where the `check_r_dependencies.yml` workflow was not properly handling branch names with `::`.
## Minor changes
- Updated the `caret` package to version 6.0-93.
- Updated the README to include information about the Open Problems team and task leaders.
- Replaced the `NuSVR` method with a faster alternative, improving performance.
## New functionality
- Added a new method for running Seuratv3 from a fork, allowing for more efficient use of resources.
- Added a new requirement to the `r_requirements.txt` file for the `bslib` package.
- Added a new requirement to the `r_requirements.txt` file for the `caret` package.
## Documentation
- Added a new section to the README to document the process of running Seuratv3 from a fork.
- Updated the README to include a list of all contributors to the Open Problems project.


# openproblems v0.5.5

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fix sampling and reindexing
- Fix docker unavailable error to include image name

## New functionality
- Require minimum celltype count for `spatial_decomposition`

## Minor changes
- Update Rcpp to 1.0.9
- Update to nf-openproblems v1.7


# openproblems v0.5.4

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed an issue where some cell types were missing from the output.


# openproblems v0.5.3

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Fixed a bug in the rctd method where cell types with fewer than 25 cells were not being used.


# openproblems v0.5.2

Note: This changelog was automatically generated from the git log.

## Bug fixes
- Handle missing function error by catching FileNotFoundError and NoSuchFunctionError instead of just RuntimeError.


# openproblems v0.5.1

Note: This changelog was automatically generated from the git log.

## Major changes
- Updated `scipy` requirement from `==1.8.*` to `>=1.8,<1.10`.
- Updated `igraph` to version `1.3.4`.

## Minor changes
- Changed the mnnpy dependency to use a patch version instead of a specific commit hash.

## Bug fixes
- Changed `docker_hash` to use the Docker API if `docker` is not available.
- Use `curl` to retrieve the Docker hash if `docker` fails.
- Fixed an issue with using `git+https` for `mnnpy`.


# openproblems v0.5.0

Note: This changelog was automatically generated from the git log.

## Minor changes

- Updated several R package dependencies.
- Updated several Python package dependencies.
- Added several new methods for spatial decomposition: RCTD, DestVI, Stereoscope.
- Added a new dataset for dimensionality reduction: Mouse hematopoietic stem cell differentiation.
- Improved documentation for tasks and datasets.

## Bug fixes

- Fixed a bug where the `lintr` package was not being installed correctly.
- Fixed a bug where the `BRANCH_PARSED` variable was not being properly sanitized in the `run_tests.yml` workflow.
- Fixed a bug in `_scanvi` and `_scvi` functions where the `max_epochs` parameter was not being passed to the `scanvi` and `scvi` functions.
- Fixed a bug in `install_renv.R` causing incorrect installation of packages from R repositories.
- Fixed an issue where the dependency upgrade script would fail to capture the output of the upgrade process.
- Fixed an issue where the dependency upgrade script would not correctly write updates to the requirements file.
- Fixed an issue where the `git_hash` function was not being called for external modules.
- Fixed a bug in `openproblems/tasks/denoising/methods/__init__.py` that prevented DCA from being used.
- Fixed a bug in `neuralee_default` where it could fail due to sparseness of data.
- Fixed a bug in `scanvi_all_genes` where the code version was not being set correctly.
- Fixed a bug in `scanvi_hvg` where the code version was not being set correctly.
- Fixed a bug in `scarches_scanvi_all_genes` where the code version was not being set correctly.
- Fixed a bug in `scarches_scanvi_hvg` where the code version was not being set correctly.

## New functionality

- Added a new denoising method called "DCA" based on a deep count autoencoder.
- Added `xgboost_log_cpm` and `xgboost_scran` methods to `openproblems.tasks.label_projection`.
- Added a new command-line interface for testing datasets, methods, and metrics.
- Added a new `install_renv.R` script to simplify the installation of the `renv` package.
- Added automated CI check to find and suggest available updates to R packages in docker images.
- Added a hash to the docker image that is based on the age of the code.
- Added data_reference to dataset metadata.
- Added `docker_hash` function to retrieve the docker image hash associated with an image.
- Added support for retrieving the docker image hash for R functions that have a defined `__r_file__`.

## Documentation

- Updated contributing guide to reflect the `main` branch as the default branch.
- Updated issue templates to reflect the `main` branch as the default branch.
- Updated pull request template to reflect the `main` branch as the default branch.


# openproblems v0.4.4

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a new docker image `openproblems-r-pytorch` for running Harmony in Python

## Major changes
- Moved `harmony` to Python-based `harmony-pytorch`

## Bug fixes
- Fixed an issue where `adata.var` was not being correctly handled in `_utils.py`
- Updated the documentation for the `openproblems-r-extras` docker image


# openproblems v0.4.3

Note: This changelog was automatically generated from the git log.

## New functionality
- Added PHATE with sqrt potential

## Bug fixes
- Fixed path to R_HOME
- Fixed Dockerfile to use R 4.2
- Minor CI fixes


# openproblems v0.4.2

Note: This changelog was automatically generated from the git log.

## Minor changes
- Run scran pooling in series, not in parallel.


# openproblems v0.4.1

Note: This changelog was automatically generated from the git log.

## New functionality
- Added `FastMNN`, `Harmony`, and `Liger` methods for batch integration.
- Added `bbknn_full_unscaled` method.
- Added Dependabot configuration for pip and GitHub Actions dependencies.

## Minor changes
- Updated dependencies: `scib`, `bbknn`, `scanorama`, `annoy`, and `mnnpy`.
- Improved the performance of several methods by pre-processing the data before running them.

## Bug fixes
- Fixed bugs in `fastMNN`, `harmony`, `liger`, `scanorama`, `scanvi`, `scvi`, `mnn`, and `combat` that caused incorrect embedding.


# openproblems v0.4.0

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a new file `workflow/generate_website_markdown.py` to generate website markdown files for all tasks and datasets.
- Updated Nextflow version to v1.5.
- Updated Nextflow version to v1.6.

## Major changes
- Added code version to the output of each method.
- Updated `nextflow` version to `v1.3`.
- Updated `nextflow` version to `v1.4`.
- Updated docker version to 20.10.15.
- Removed Docker setup from CI workflow.
- Updated Python version to 3.8.13.

## Minor changes
- Updated dependencies for the Docker images.
- Updated pre-commit hooks to include `requirements-txt-fixer`.
- Updated Nextflow workflow to version 1.4.
- Updated the location of method versions in the results directory.
- Updated the Tower action ID.

## Bug fixes
- Fixed a bug where Docker images were not properly pushed to Docker Hub.
- Updated `requirements.txt` files to fix dependency conflicts.
- Removed unnecessary dependencies from CI workflows to reduce disk space usage on GitHub runners.


# openproblems v0.3.5

Note: This changelog was automatically generated from the git log.

## New functionality
- Added new integration methods: BBKNN, Combat, FastMNN feature, FastMNN embed, Harmony, Liger, MNN, Scanorama feature, Scanorama embed, Scanvi, Scvi
- Added new metrics: graph_connectivity, iso_label_f1, nmi
- Added _utils.py with functions: hvg_batch, scale_batch
- Added `run_bbknn` function.
- Added a test for the trustworthiness metric, which now passes for sparse matrices.
- Added a test for the density preservation metric, which now passes against densmap for a reasonable degree of similarity.
- Added tests for all methods and metrics.
- Added a new workflow to automatically delete untagged images from the OpenProblems ECR repository.
- Added a new workflow to process results and create a PR to update the OpenProblems benchmark.
- Added support for running tests with the `process` extra in `setup.py`.
- Added `densmap` dimensionality reduction method.
- Added `neuralee` dimensionality reduction method.
- Added `alra` denoising method.
- Added `scarches_scanvi` label projection method.
- Added `bbknn` batch integration graph method.
- Added `beta` regulatory effect prediction method.
- Added a new `invite-contributors.yml` file to the repository.

## Major changes
- The `test_methods.py` file has been simplified by removing unused arguments.
- The `test_metrics.py` file has been simplified by removing unused arguments.
- The `test_utils/docker.py` file has been modified to allow specifying the docker image as a decorator argument.
- Updated Nextflow version to 22.04.0.
- Modified the processing of Nextflow results to save them in a temporary directory.
- Modified `workflow/parse_nextflow.py` to parse results from Nextflow runs.
- Modified `.github/workflows/run_tests.yml` to cancel previous runs when a new commit is pushed.

## Minor changes
- Removed `.nextflow`, `scratch/`, `openproblems/results/` and `openproblems/work/` from `.gitignore`.
- Updated `CONTRIBUTING.md`
- Methods should not edit `adata.obsm["train"]` or `adata.obsm["test"]`.
- Redirects stdout to stderr when running subcommands to ensure that output is printed correctly.
- Updated CI workflow to skip running tests on push if they failed on the `run_tester` job, unless the branch name starts with `test_benchmark`.
- Refactored the Neuralee method to use a separate function for embedding.
- Improved performance by using a default value for `maxit` in `fine_tune_kwargs`.
- Removed unnecessary code for storing raw counts in the `neuralee_default` method.

# openproblems v0.3.4

Note: This changelog was automatically generated from the git log.

## New functionality
- Added CeNGEN, Tabula Muris Senis, and Pancreas datasets to the label_projection task.
- Added scANVI and scArches+scANVI methods to the label_projection task.
- Added majority_vote and random_labels baseline methods to the label_projection task.
- Added new methods: densMAP, NeuralEE, scvis
- Added new metrics: NN Ranking (continuity, co-KNN size, co-KNN AUC, Local continuity meta criterion, Local property metric, Global property metric)
- Added pre-processing function: log_cpm_hvg()
- Added support for custom pre-processing functions
- Added support for variants of methods
- Added a new batch integration task.
- Added a batch integration graph subtask.
- Added a batch integration embedding subtask.
- Added a batch integration corrected feature matrix subtask.
- Added ivis method for dimensionality reduction to openproblems.
- Added self-hosted runner support for `run_benchmark` workflow using Cirun.io
- Added a `--test` flag to the `run` subcommand, allowing for running a test version of a method.
- Added `test_load_dataset` to `test/test__load_data.py` to test loading and caching of datasets.
- Added `test_method` to `test/test_methods.py` to test application of methods.
- Added `test_trustworthiness_sparse` to `test/test_metrics.py` to test trustworthiness metric on sparse data.
- Added `test_density_preservation_matches_densmap` to `test/test_metrics.py` to test density preservation metric against densmap.
- Updated `test/utils/docker.py` to allow specifying the docker image as the last argument.
- Added `--test` flag to `run` subcommand to run the test version of a method.
- Added Docker image building to `run_tests.yml`.
- Added a new workflow to process Nextflow results
- Added a new workflow to run tests and benchmarks
- Added support for running benchmarks from tags
- Added support for running benchmarks from forks
- Added `openproblems-cli` command to run test-hash


# openproblems v0.3.3

Note: This changelog was automatically generated from the git log.

## New functionality
- Added support for balanced SCOT alignment.

## Minor changes
- Updated the workflow to store benchmark results in `/tmp`.

## Bug fixes
- Fixed the parsing and committing of benchmark results on tag.
- Fixed the Github Actions badge link.
- Fixed the coverage badge.
- Fixed the benchmark commit.
- Ignored AWS warning and cleaned up S3 properly.
- Updated the workflow to continue on error for forks.


# openproblems v0.3.2

Note: This changelog was automatically generated from the git log.

## New functionality
- Added trustworthiness metric to the dimensionality reduction task.
- Added density preservation metric.
- Added several metrics based on nearest neighbor ranking: continuity, co-KNN size, co-KNN AUC, local continuity meta criterion, local property, global property.
- Added mouse blood data from Olsson et al. (2016) Nature to the `openproblems` dataset collection.
- Added a test mode to the `load_olsson_2016_mouse_blood` function.
- Added a dataset function for the `mouse_blood_olssen_labelled` dataset in the `openproblems.tasks.dimensionality_reduction.datasets` module.
- Added ALRA denoising method.
- Added support for the Single Cell Optimal Transport (SCOT) method for multimodal data integration.
- SCOT implements Gromov-Wasserstein optimal transport to align single-cell multi-omics data.
- Added four variations of SCOT:
+    - sqrt CPM unbalanced
+    - sqrt CPM balanced
+    - log scran unbalanced
+    - log scran balanced
- Each variation implements different normalization strategies for the input data.
- Added `scot` method to `openproblems.tasks.multimodal_data_integration.methods`.
- Added pre-processing to the `dimensionality_reduction` task.
- Added pre-processing to all `dimensionality_reduction` methods.
- Added Wagner_2018_zebrafish_embryo_CRISPR dataset loader
- Added PR review checklist to the pull request template.
- Added `cmake==3.18.4` to the `docker/openproblems-python-extras/requirements.txt` file.
- Added `--version` flag to print the version.
- Added `--test-hash` flag to print the current hash.
- Added basic help message.
- Added `install_renv.R` script for installing R packages in Docker images.
- Added `docker/.version` file to track Docker image version.
- Added a new docker image for running GitHub Actions.
- Added a new utils.git module to determine which tasks have changed relative to base/main.
- Added support for running benchmark tests on tags.
- Added a test directory for use in the workflow.


# openproblems v0.3.1

Note: This changelog was automatically generated from the git log.

## New functionality
- Added chromatin potential task
- Added PHATE to the dimensional_reduction task.
- Added support for testing docker builds on a separate branch.
- Added support for building images and pushing them to docker hub.
- Added support for writing methods in R using `scprep`'s `RFunction` class.
- Added a CLI interface to `openproblems`.
- Added `f1_micro` metric.
- Added `mlp_log_cpm` and `mlp_scran` methods for label projection.
- Added `pancreas_batch` and `pancreas_random` datasets for label projection.
- Added `f1` metric for label projection.
- Added metadata to methods and metrics.
- Added `openproblems.tools.decorators` for decorating methods and metrics.
- Added `openproblems.tools.normalize` for common normalization functions.
- Added methods for `logistic_regression`, `mlp`, `harmonic_alignment`, `mnn`, and `procrustes`.
- Added metrics for `accuracy`, `f1`, `knn_auc`, and `mse`.
- Added `openproblems.version` to provide package version.
- Added `dataset` decorator for registering datasets.
- Added `tools.decorators.profile` decorator to measure memory usage and runtime of methods.
- Added `tools.normalize` module to provide normalization functions.
- Added `tools.decorators.normalizer` decorator to normalize data prior to applying methods.
- Added a new "data loader" component that loads data in a way that's formatted correctly for a given task.
- Added CITE-seq Cord Blood Mononuclear Cells dataset.
- Added snakemake support for automatic evaluation.
- Added zebrafish data to label projection task.
- Added a new task, "Link gene expression with chromatin accessibility"
- Added a new dataset, "sciCAR Mouse Kidney with cell clusters"
- Added a new method, "BETA"
- Added a new metric, "Correlation between RNA and ATAC"
- Added a new task, "Dimensional reduction"
- Added human blood dataset from Nestorowa et al. Blood. 2016
- Added 10x PBMC dataset
- Added `load_10x_5k_pbmc` function to load the 10x 5k PBMC dataset.


# openproblems v0.2.1

Note: This changelog was automatically generated from the git log.

## New functionality
- Added MLP method for label projection task.
- Added pancreas data loading to label projection task.

## Minor changes
- Updated black.
- Updated test version of pancreas_batch to have test data.
- Added random pancreas train data.

## Bug fixes
- Fixed zebrafish code duplication.
- Fixed pancreas import location.
- Fixed bug in zebrafish data.
- Fixed bug in pancreas import.
- Removed normalization from loader.
- Removed dummy and cheat metrics/datasets.
- Removed excess covariates from pancreas dataset.


# openproblems v0.2

Note: This changelog was automatically generated from the git log.

## New functionality
- Added zebrafish label projection task

## Major changes
- Moved scIB, rpy2, harmonicalignment, and mnnpy to optional dependencies

## Minor changes
- Improved n_components fix
- Moved URL into function for neater namespace

## Bug fixes
- Fixed n_svd for truncatedSVD
- Fixed data loader
- Fixed n_pca problem
- Scaled without mean if sparse
- Scaled data for regression
- Added check to ensure that data has nonzero size


# openproblems v0.1

Note: This changelog was automatically generated from the git log.

## New functionality
- Added a results page to the website.
- Added a new zebrafish dataset to the openproblems library.
- Added netlify.toml to deploy website.

## Documentation
- Updated documentation to reflect new features and datasets.

## Major changes
- Bumped version to 0.1.

## Minor changes
- Improved the website's home menu link.
- Improved website links.
- Updated website's hero and social links.
- Updated website's task cards.
- Updated the website's demo.
- Improved website's frontmatter.
- Separated frontmatter from content in website's Markdown files.
- Fixed black syntax.
- Excluded website from black.
- Updated website content to display results.
- Updated the Travis CI configuration to exclude website from black.

## Bug fixes
- Fixed zebrafish data loader.


# openproblems v0.0.3

Note: This changelog was automatically generated from the git log.

## New functionality
- Added harmonic alignment method.
- Added scicar datasets.
- Added logistic regression methods.
- Added ability to normalize obsm.
- Added test suite.
- Added normalization tools.

## Documentation
- Updated documentation to reflect normalization changes.

## Major changes
- Migrated normalizations to openproblems.tools.normalize.
- Updated dataset specification to require normalization in methods.
- Removed zebrafish dataset.
- Moved dataset test spec.
- Removed "mode2_raw" and "raw" from datasets.
- Added test dataset spec.
- Consolidated scicar datasets.
- Migrated references to github repo.

## Minor changes
- Improved sparse array equality test.
- Improved sparse inequality check.
- Increased test data size.
- Normalized mode2.
- Fixed decorator.
- Used uns.
- Used functools.wraps.
- Updated name of log_scran_pooling function.
- Fixed storing normalization results.
- Fixed zebrafish load caching.
- Fixed zebrafish test.
- Added normalization functions.
- Updated logistic regression function to work with anndata properly.
- Fixed cheat method.
- Fixed git upload.
- Fixed Travis CI.
- Fixed harmonic alignment import.
- Increased test coverage.

## Bug fixes
- Bugfix harmonic_alignment, closes #4.
- Bugfix harmonic alignment import.
- Normalized data inside methods, closes #19.
- Fix storing normalization results.
- Fixed zebrafish test.
- Fix zebrafish load caching.
- Fix decorator.
- Fix cheat method.
- Don't check for raw data -- we are no longer normalizing.


# openproblems v0.0.2

Note: This changelog was automatically generated from the git log.

## New functionality

- Added dummy dataset to `openproblems/data`
- Added `load_dummy` function to `openproblems/data`
- Added `loader` decorator to `openproblems/data`
- Added loading functions for sciCAR datasets to `openproblems/data/scicar`
- Added `scicar_cell_lines` dataset to `openproblems/tasks/multimodal_data_integration/datasets`
- Added `scicar_mouse_kidney` dataset to `openproblems/tasks/multimodal_data_integration/datasets`
- Added `dummy` dataset to `openproblems/tasks/label_projection/datasets`

## Major changes

- Changed data structure for multimodal data integration tasks in `openproblems/tasks/multimodal_data_integration`
- Bumped version to 0.0.2 in `openproblems/version.py`
- Modified the way to run `evaluate.sh` in `.travis.yml`
- Added `chmod +x evaluate.sh` to `.travis.yml`

## Documentation

- Added documentation for adding a dataset to a task in `README.md`
- Added documentation for dataset loading in `README.md`
- Added documentation for adding a new dataset in `README.md`
- Updated documentation in `openproblems/tasks/multimodal_data_integration/README.md`
- Updated documentation in `openproblems/version.py`

# openproblems v0.0.1

First release of OpenProblems.

methods, 1 metric)
* Multimodal data integration (2 datasets, 2 methods, 2 metrics)
