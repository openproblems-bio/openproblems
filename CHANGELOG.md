
# openproblems-v2 0.1.0

## common

### NEW FUNCTIONALITY

* `extract_scores`: Summarise a metrics output tsv.

* `dataset_concatenate`: Concatenate N AnnData datasets.

* Created test data `resources_test/pancreas` with `src/common/resources_test_scripts/pancreas.sh`.

* `list_git_shas`: create list of latest commit hashes of all files in repo.

* `get_api_info`: extract api info from tasks

* `get_method_info`: extract method info from config yaml

* `get_metric_info`: extract metric info from config yaml

* `check_migration_status`: compare git shas from methods with v1

* `get_results`: extract benchmark scores 

* `get_task_info`: extract task info


## label_projection

### NEW FUNCTIONALITY

* `api/anndata_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, method and metric components.

* `split_dataset`: Added a component for splitting raw datasets into task-ready dataset objects.

* `resources_test/label_projection/pancreas` with `src/label_projection/resources_test_scripts/pancreas.sh`.

### V1 MIGRATION

* `methods/knn`: Migrated from v1.

* `methods/logistic_regression`: Migrated from v1.

* `methods/mlp`: Migrated from v1.

* `methods/scanvi`: Migrated and adapted from v1.

* `methods/seurat_transferdata`: Migrated and adapted from v1.

* `methods/xgboost`: Migrated from v1.

* `control_methods/majority_vote`: Migrated from v1.

* `control_methods/random_labels`: Migrated from v1.

* `control_methods/true_labels`: Migrated from v1.

* `metric/accuracy`: Migrated from v1.

* `metric/f1`: Migrated from v1.

## datasets

### NEW FUNCTIONALITY

* `workflows/process_openproblems_v1`: Fetch and process legacy OpenProblems v1 datasets

* `normalization/log_cpm`: A log CPM normalization method.

* `normalization/log_scran_pooling`: A log scran pooling normalization method.

* `normalization/sqrt_cpm`: A sqrt CPM normalization method.

* `normalization/l1_sqrt`: A scaled L1 sqrt normalization. extracted from Alra method in the denoising task from v1

* `subsample`: Subsample an h5ad file.

### V1 MIGRATION

* `loaders/openproblems_v1`: Fetch a dataset from OpenProblems v1

## denoising

### NEW FUNCTIONALITY

* `api/anndata_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, method and metric components.

* `split_dataset`: Added a component for splitting raw datasets into task-ready dataset objects.

* `resources_test/denoising/pancreas` with `src/denoising/resources_test_scripts/pancreas.sh`.

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

* split_dataset also removes unnecessary data from train and test datasets not needed by the methods and metrics.

## Dimensionality reduction

### New functionality
* `api/anndata_*`: Created a file format specifications for the h5ad files throughout the pipeline.

* `api/comp_*`: Created an api definition for the split, control method, method and metric components.

* `split_dataset`: Added a component for splitting raw datasets into task-ready dataset objects.

* `control_methods`: Added a component for baseline methods specifically.

* `resources_test/dimensionality_reduction/pancreas` with `src/dimensionality_reduction/resources_test_scripts/pancreas.sh`.

* Added `variant` key to config files to store variants (different input parameters) of every component.

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

* `metrics/rmse`: Migrated from v1, but will likely be removed.

* `metrics/trustworthiness`: Migrated from v1, but will likely be removed.

* `metrics/density_preservation`: Migrated from v1.

* `metrics/coranking`: Migrated from v1. This script originally called `nn_ranking.py` and written in Python.

### Changes from V1

* Raw counts and normalized expression data is stored in `.layers["counts"]` and `.layers["normalized"]`, respectively,
  instead of in `.X`.
  
* A `split_dataset` has been implemented to make a distinction between the data a method is allowed to see
  (here called the train data) and what a metric is allowed to see (here called the test data).

* `methods/ivis` had originally been removed from the v1 (temporarily) but has been added back to the v2.

* The metrics as defined in the `nn_ranking.py` script have been documented and refactored into an R
  component `metrics/coranking`.

* `metrics/rmse` should be removed because RMSE metrics don't really make sense here.

* `metrics/trustworthiness` should be removed because it is already included in `metrics/coranking`.
