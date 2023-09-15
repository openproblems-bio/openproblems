
# openproblems-v2 0.1.0

## general

### NEW FUNCTIONALITY

* Updated all current tasks in v2 to latest changes in OP v1 (PR #214)

### MAJOR CHANGES

* Relocate task directories to new `src/tasks/` location (PR #142).

* Update Docker images to our base images; `ghcr.io/openproblems-bio/base-python`
  and `ghcr.io/openproblems-bio/base-r` (PR #168).

* Update batch integration docker images to OpenProblems base images (PR #171).
  
* Changed default normalization CPM to CP10k (PR #214)

### MINOR CHANGES

* Update test scripts (PR #143).

* Update "baseline" to "control" (PR #146).

* Add task image thumbnails to api (PR #231).

### BUG FIXES

* `dimensionality_reduction/methods/tsne`: Use GitHub version of MulticoreTSNE.

* `label_projection/methods/seurat_transferdata`: Temporarily disable component as it appears to not be working (PR #206).

* Remove the ns-list action for workflows in integration test (PR #208)


## common

### NEW FUNCTIONALITY

* `extract_scores`: Summarise a metrics output tsv.

* Created test data `resources_test/pancreas` with `src/common/resources_test_scripts/pancreas.sh`.

* `get_api_info`: Extract api info from tasks.

* `get_method_info`: Extract method info from config yaml.

* `get_metric_info`: Extract metric info from config yaml.

* `get_results`: Extract benchmark scores.

* `get_task_info`: Extract task info.

* `comp_tests`: Common unit tests that can be used by all tasks.

* `check_dataset_schema`: Check if the dataset used has the required fields defined in the api `file_*.yaml` files.
  
* `Create_component`: Creates a template folder with a viash config and script file depending on the task api.

### MINOR CHANGES

* Refactor and standardize metric and method info fields (PR #99).

* Add url check to method and metric unit test (PR #160).

* Add library.bib file check to component unit test (PR #167)

### BUG FIXES

* fix typos in metric and common defenition schemas (PR #212)

## migration

### NEW FUNCTIONALITY

* `list_git_shas`: create list of latest commit hashes of all files in repo.

* `check_migration_status`: compare git shas from methods with v1

## datasets

### NEW FUNCTIONALITY

* `workflows/process_openproblems_v1`: Fetch and process legacy OpenProblems v1 datasets, whilst adding extra information to the `.uns`.

* `normalization/log_cpm`: A log CPM normalization method.

* `normalization/log_scran_pooling`: A log scran pooling normalization method.

* `normalization/sqrt_cpm`: A sqrt CPM normalization method.

* `normalization/l1_sqrt`: A scaled L1 sqrt normalization. extracted from Alra method in the denoising task from v1

* `subsample`: Subsample an h5ad file. Allows keeping observations from specific batches and celltypes, 
  also allows keeping certain features.

* `resources_test_scripts`: Scripts to create test_resources for local development with "pancreas", "pancreas_tasks" and "multimodal".

### V1 MIGRATION

* `loaders/openproblems_v1`: Fetch a dataset from OpenProblems v1, whilst adding extra information to the `.uns`.

* `loaders/openproblems_v1_multimodal`: Fetch a multimodal dataset from OpenProblems v1, whilst adding extra information to the `.uns`.

## batch_integration

### NEW FUNCTIONALITY

* `api/file_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the process, method and metric components.

* `process_dataset`: Added a component for processing common datasets into task-ready dataset objects.

* `resources_test/label_projection/pancreas` with `src/tasks/label_projection/resources_test_scripts/pancreas.sh`.

* `workflows/run`: Added nf-tower test script (PR #205).

* `metrics/lisi`: Added a component for cLISI and iLISI graph metrics from scib (PR #213).

### V1 MIGRATION

* Removed the separate subtask specific subfolders. The "subtask" is added to the config.

* `control_methods/no_integration_batch`: Migrated from v1 embedding.

* `control_methods/random_embed_cell`: Migrated from v1 embedding.

* `control_methods/random_embed_cel_jitter`: Migrated from v1 embedding.

* `control_methods/random_integration`: Migrated from v1 graph.

* `methods/bbknn`: Migrated from v1 graph.

* `methods/combat`: Migrated from v1 feature.

* `methods/scanorama_embed`: Migrated from v1 embedding.

* `methods/scanorama_feature`: Migrated from v1 feature.

* `methods/scvi`: Migrated from v1 embedding.

* `metrics/asw_batch`: Migrated from v1 embedding.

* `metrics/asw_label`: Migrated from v1 embedding.

* `metrics/cell_cycle_conservation`: Migrated from v1 embedding.

* `metrics/clustering_overlap`: Migrated from v1 graph NMI & ARI.

* `metrics/graph_connectivity`: Migrated from v1 graph.

* `metrics/hvg_overlap`: Migrated from v1 feature.

* `metrics/isolated_label_asw`: Migrated from v1 embedding.

* `metrics/isolated_label_f1`: Migrated from v1 graph.

* `metrics/kbet`: Migrated from v1 embedding.

* `metrics/pcr`: Migrated from v1 embedding.

### MINOR CHANGES

* Removed the `.uns["output_type"]` field from output anndata in methods and control methods. (PR #205)

## label_projection

### NEW FUNCTIONALITY

* `api/file_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the process, method and metric components.

* `process_dataset`: Added a component for processing common datasets into task-ready dataset objects.

* `resources_test/label_projection/pancreas` with `src/tasks/label_projection/resources_test_scripts/pancreas.sh`.

* * `workflows/run`: Added nf-tower test script. (PR #205)

### V1 MIGRATION

* `methods/knn`: Migrated from v1.

* `methods/logistic_regression`: Migrated from v1.

* `methods/mlp`: Migrated from v1.

* `methods/scanvi`: Migrated and adapted from v1.

* `methods/scanvi_scarches`: Migrated and adapted from v1.

* `methods/seurat_transferdata`: Migrated and adapted from v1.

* `methods/xgboost`: Migrated from v1.

* `control_methods/majority_vote`: Migrated from v1.

* `control_methods/random_labels`: Migrated from v1.

* `control_methods/true_labels`: Migrated from v1.

* `metric/accuracy`: Migrated from v1.

* `metric/f1`: Migrated from v1.

## denoising

### NEW FUNCTIONALITY

* `api/file_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, method and metric components.

* `process_dataset`: Added a component for processing common datasets into task-ready dataset objects.

* `resources_test/denoising/pancreas` with `src/tasks/denoising/resources_test_scripts/pancreas.sh`.
  
* `workflows/run`: Added nf-tower test script. (PR #205)

### V1 MIGRATION

* `control_methods/no_denoising`: Migrated from v1. Extracted from baseline method

* `control_methods/perfect_denoising`: Migrated from v1.Extracted from baseline method

* `methods/alra`: Migrated from v1. Changed from python to R and uses lg_cpm normalised data instead of L1 sqrt

* `methods/dca`: Migrated and adapted from v1.

* `methods/knn_smoothing`: Migrated and adapted from v1.

* `methods/magic`: Migrated from v1.

* `metrics/mse`: Migrated from v1.

* `metrics/poisson`: Migrated from v1.

### Changes from V1

* Anndata layers are used to store data instead of obsm
  
* extended the use of sparse data in methods unless it was not possible

* process_dataset also removes unnecessary data from train and test datasets not needed by the methods and metrics.

## Dimensionality reduction

### New functionality
* `api/file_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, control method, method and metric components.

* `process_dataset`: Added a component for processing common datasets into task-ready dataset objects.

* `control_methods`: Added a component for baseline methods specifically.

* `resources_test/dimensionality_reduction/pancreas` with `src/tasks/dimensionality_reduction/resources_test_scripts/pancreas.sh`.

* Added `variant` key to config files to store variants (different input parameters) of every component.
  
* `workflows/run`: Added nf-tower test script. (PR #205)

### V1 migration
* `control_methods/true_features`: Migrated from v1. Extracted from baseline method `True Features`.

* `control_methods/random_features`: Migrated from v1. Extracted from baseline method `Random Features`.

* `methods/umap`: Migrated from v1.

* `methods/ivis`: Migrated from v1.

* `methods/tsne`: Migrated and adapted from v1.

* `methods/densmap`: Migrated and adapted from v1.

* `methods/phate`: Migrated from v1.

* `methods/pca`: Migrated from v1.

* `methods/neuralee`: Migrated from v1.

* `metrics/distance_correlation`: Migrated from v1, but will likely be removed.

* `metrics/trustworthiness`: Migrated from v1, but will likely be removed.

* `metrics/density_preservation`: Migrated from v1.

* `metrics/coranking`: Migrated from v1. This script originally called `nn_ranking.py` and written in Python.

### Changes from V1

* Raw counts and normalized expression data is stored in `.layers["counts"]` and `.layers["normalized"]`, respectively,
  instead of in `.X`.
  
* A `process_dataset` has been implemented to make a distinction between the data a method is allowed to see
  (here called the train data) and what a metric is allowed to see (here called the test data).

* `methods/ivis` had originally been removed from the v1 (temporarily) but has been added back to the v2.

* The metrics as defined in the `nn_ranking.py` script have been documented and refactored into an R
  component `metrics/coranking`.

* `metrics/rmse` should be removed because RMSE metrics don't really make sense here.

* `metrics/trustworthiness` should be removed because it is already included in `metrics/coranking`.


## match_modalities (PR #201)

### New functionality

* `api/file_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, control method, method and metric components.

* `process_dataset`: Added a component for processing common datasets into task-ready dataset objects.

* `control_methods`: Added a component for baseline methods specifically.

* `resources_test/dimensionality_reduction/pancreas` with `src/tasks/dimensionality_reduction/resources_test_scripts/pancreas.sh`.

* Added `variant` key to config files to store variants (different input parameters) of every component.
  
* `workflows/run`: Added nf-tower test script.

### V1 migration

* `control_methods/true_features`: Migrated from v1. Extracted from baseline method `True Features`.

* `control_methods/random_features`: Migrated from v1. Extracted from baseline method `Random Features`.

* `methods/harmonic_alignment`: Migrated from v1.

* `methods/mnn`: Migrated from v1.

* `methods/procrustes`: Migrated from v1.

* `metrics/knn_auc`: Migrated from v1.

* `metrics/mse`: Migrated from v1.


### Changes from V1

* `methods/scot`: Add new scot method.

* Raw counts and normalized expression data is stored in `.layers["counts"]` and `.layers["normalized"]`, respectively,
  instead of in `.X`.

* The methods and metrics now take 2 modal datasets as input instead of 1.


