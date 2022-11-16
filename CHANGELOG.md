
# openproblems-v2 0.1.0

## common

### NEW FUNCTIONALITY

* `extract_scores`: Summarise a metrics output tsv.

* `dataset_concatenate`: Concatenate N AnnData datasets.

* Created test data `resources_test/pancreas` with `src/common/resources_test_scripts/pancreas.sh`.


## label_projection

### NEW FUNCTIONALITY

* API: Created an explicit api definition for the split, method and metric components.

* `data_processing/split`: Added a component for splitting raw datasets into task-ready dataset objects.

* `resources_test/label_projection/pancreas` with `src/label_projection/resources_test_scripts/pancreas.sh`.

### V1 MIGRATION

* `methods/knn_classifier`: Migrated from v1.

* `methods/logistic_regression`: Migrated from v1.

* `methods/mlp`: Migrated from v1.

* Temporarily disable `scanvi` / `scarches_scanvi`.

* `metric/accuracy`: Migrated from v1.

* `metric/f1`: Migrated from v1.

## datasets

### NEW FUNCTIONALITY

* `workflows/process_openproblems_v1`: Fetch and process legacy OpenProblems v1 datasets

* `normalization/log_cpm`: A log CPM normalization method.

* `normalization/log_scran_pooling`: A log scran pooling normalization method.

* `normalization/sqrt_cpm`: A sqrt CPM normalization method.

* `subsample`: Subsample an h5ad file.

### V1 MIGRATION

* `loaders/openproblems_v1`: Fetch a dataset from OpenProblems v1
