# openproblems-v2 0.1.0

## common

### NEW FUNCTIONALITY

* `dataset_loader/download`: Download an AnnData dataset from a URL.

* `extract_scores`: Summarise a metrics output tsv.

* `dataset_concatenate`: Concatenate N AnnData datasets.

* `subsample`: Subsample an anndata file.

* Created test data `resources_test/pancreas` with `src/common/resources_test_scripts/pancreas.sh`.

* Created normalization method `common/normalise_log_cpm`.

* Created normalization method `common/normalise_log_scran_pooling`.

## label_projection

### NEW FUNCTIONALITY

* API: Created an explicit api definition for the censor, method and metric components.

* Created censoring component `data_processing/censoring`.

* Created test data `resources_test/label_projection/pancreas` with `src/label_projection/resources_test_scripts/pancreas.sh`.

### V1 MIGRATION

* Ported method `methods/knn_classifier`.

* Ported method `methods/logistic_regression`.

* Ported method `methods/mlp`.