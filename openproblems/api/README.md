# openproblems-cli

This CLI provides an interface to the `openproblems` Python package from the command line, primarily for the testing and evaluation workflow.

## Usage

```
usage: openproblems-cli [-h] [--parallel]
                        {tasks,list,image,load,run,evaluate} ...

Open Problems for Single Cell Analysis command-line interface

positional arguments:
  {tasks,list,image,load,run,evaluate}
    tasks               List tasks
    list                List datasets, methods and metrics
    image               Fetch a Docker image associated with a function
    load                Load a dataset for a task
    run                 Run a method on a dataset
    evaluate            Evaluate a metric on a method

optional arguments:
  -h, --help            show this help message and exit
  --parallel, -p        Run tasks in parallel. This prevents deletion of the
                        cache
```

## Example (without docker)
Running the CLI requires commands to be run in a specific order: `load` -> `run` -> `evaluate`.

For example:
```
# Download a task-specific dataset and save it to `dataset.h5ad`
openproblems-cli load --task label_projection --output dataset.h5ad pancreas_batch
# Run a method on a datasets and save output to `method.h5ad`
openproblems-cli run --task label_projection --input dataset.h5ad --output method.h5ad logistic_regression_log_cpm
# Evaluate the performance of a previously run method using the `accuracy` metric
openproblems-cli evaluate --task label_projection --input method.h5ad accuracy
```

You can list available tasks using `openproblems-cli tasks`
```
> openproblems-cli tasks
denoising
dimensionality_reduction
label_projection
multimodal_data_integration
regulatory_effect_prediction
```

You can then list the avaiable datasets, methods, and metrics for a partiular task using `openproblems-cli list --[datasets|methods|metrics] --task [task_name]`
```
> openproblems-cli list --datasets --task label_projection
pancreas_batch
pancreas_random
zebrafish_labels
zebrafish_random

> openproblems-cli list --methods --task label_projection
knn_classifier_log_cpm
knn_classifier_scran
logistic_regression_log_cpm
logistic_regression_scran
mlp_log_cpm
mlp_scran

> openproblems-cli list --metrics --task label_projection
accuracy
f1
f1_micro
```

The output of these commands are allowable arguments to the respective `load`, `run`, and `evaluate` commands.

### Sample output

```
$ openproblems-cli tasks
chromatin_potential
denoising
dimensional_reduction
label_projection
multimodal_data_integration
$ openproblems-cli list --datasets --task label_projection
pancreas_batch
pancreas_random
zebrafish_labels
zebrafish_random
$ openproblems-cli load --task label_projection --output dataset.h5ad pancreas_batch
$ openproblems-cli list --methods --task label_projection
logistic_regression_log_cpm
logistic_regression_scran
mlp_log_cpm
mlp_scran
$ openproblems-cli run --task label_projection --input dataset.h5ad --output method.h5ad logistic_regression_log_cpm
$ openproblems-cli list --metrics --task label_projection
$ openproblems-cli evaluate --task label_projection --input method.h5ad accuracy
0.9521233432512848
```

## Example (with docker)

```
openproblems-cli tasks
openproblems-cli list --datasets --task label_projection
openproblems-cli image --datasets --task label_projection pancreas_batch
docker run -dt openproblems-cli load --task label_projection --output dataset.h5ad pancreas_batch
openproblems-cli list --methods --task label_projection
openproblems-cli image --methods --task label_projection logistic_regression_scran
openproblems-cli run --task label_projection --input dataset.h5ad --output method.h5ad logistic_regression_log_cpm
openproblems-cli list --metrics --task label_projection
openproblems-cli image --metrics --task label_projection accuracy
openproblems-cli evaluate --task label_projection --input method.h5ad accuracy
```

### Sample output

```
$ openproblems-cli tasks
chromatin_potential
denoising
dimensional_reduction
label_projection
multimodal_data_integration
$ openproblems-cli list --datasets --task label_projection
pancreas_batch
pancreas_random
zebrafish_labels
zebrafish_random
$ openproblems-cli image --datasets --task label_projection pancreas_batch
openproblems
$ docker run -dt singlecellopenproblems/openproblems openproblems-cli load --task label_projection --output dataset.h5ad pancreas_batch
$ openproblems-cli list --methods --task label_projection
logistic_regression_log_cpm
logistic_regression_scran
mlp_log_cpm
mlp_scran
$ openproblems-cli image --methods --task label_projection logistic_regression_scran
openproblems-r-base
$ docker run -dt singlecellopenproblems/openproblems-r-base openproblems-cli run --task label_projection --input dataset.h5ad --output method.h5ad logistic_regression_log_cpm
$ openproblems-cli list --metrics --task label_projection
accuracy
f1
f1_micro
$ openproblems-cli image --metrics --task label_projection accuracy
openproblems
$ docker run -dt singlecellopenproblems/openproblems openproblems-cli evaluate --task label_projection --input method.h5ad accuracy
0.9521233432512848
```
