Contributing to OpenProblems
================

- <a href="#code-of-conduct" id="toc-code-of-conduct">Code of conduct</a>
- <a href="#requirements" id="toc-requirements">Requirements</a>
- <a href="#quick-start" id="toc-quick-start">Quick start</a>
- <a href="#project-structure" id="toc-project-structure">Project
  structure</a>
- <a href="#adding-a-viash-component"
  id="toc-adding-a-viash-component">Adding a Viash component</a>
- <a href="#running-a-component-from-cli"
  id="toc-running-a-component-from-cli">Running a component from CLI</a>
- <a href="#building-a-component" id="toc-building-a-component">Building a
  component</a>
- <a href="#unit-testing-a-component"
  id="toc-unit-testing-a-component">Unit testing a component</a>
- <a href="#more-information" id="toc-more-information">More
  information</a>

[OpenProblems](https://openproblems.bio) is a community effort, and
everyone is welcome to contribute. This project is hosted on
[github.com/openproblems-bio/openproblems-v2](https://github.com/openproblems-bio/openproblems-v2).

## Code of conduct

We as members, contributors, and leaders pledge to make participation in
our community a harassment-free experience for everyone, regardless of
age, body size, visible or invisible disability, ethnicity, sex
characteristics, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance,
race, caste, color, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open,
welcoming, diverse, inclusive, and healthy community.

Our full [Code of Conduct](CODE_OF_CONDUCT.md) is adapted from the
[Contributor Covenant](https://www.contributor-covenant.org), version
2.1.

## Requirements

To use this repository, please install the following dependencies:

- Bash
- Java (Java 11 or higher)
- Docker (Instructions [here](https://docs.docker.com/get-docker/))
- Nextflow (Optional, though [very easy to
  install](https://www.nextflow.io/index.html#GetStarted))

## Quick start

The `src/` folder contains modular software components for running a
modality alignment benchmark. Running the full pipeline is quite easy.

**Step 0, fetch Viash and Nextflow**

``` bash
mkdir $HOME/bin
curl -fsSL get.viash.io | bash -s -- --bin $HOME/bin --tools false
curl -s https://get.nextflow.io | bash; mv nextflow $HOME/bin
```

Make sure that Viash and Nextflow are on the \$PATH by checking whether
the following commands work:

``` bash
viash -v
nextflow -v
```

    viash 0.6.6 (c) 2020 Data Intuitive
    nextflow version 22.10.4.5836

**Step 1, download test resources:** by running the following command.

``` bash
viash run src/common/sync_test_resources/config.vsh.yaml
```

    Completed 256.0 KiB/7.2 MiB (302.6 KiB/s) with 6 file(s) remaining
    Completed 512.0 KiB/7.2 MiB (595.8 KiB/s) with 6 file(s) remaining
    Completed 768.0 KiB/7.2 MiB (880.3 KiB/s) with 6 file(s) remaining
    Completed 1.0 MiB/7.2 MiB (1.1 MiB/s) with 6 file(s) remaining    
    Completed 1.2 MiB/7.2 MiB (1.3 MiB/s) with 6 file(s) remaining
    ...

**Step 2, build all the components:** in the `src/` folder as standalone
executables in the `target/` folder. Use the `-q 'xxx'` parameter to
build a subset of components in the repository.

``` bash
viash ns build --query 'label_projection|common' --parallel --setup cachedbuild
```

    In development mode with 'dev'.
    Exporting split_dataset (label_projection) =docker=> target/docker/label_projection/split_dataset
    Exporting accuracy (label_projection/metrics) =docker=> target/docker/label_projection/metrics/accuracy
    Exporting random_labels (label_projection/control_methods) =docker=> target/docker/label_projection/control_methods/random_labels
    [notice] Building container 'label_projection/control_methods_random_labels:dev' with Dockerfile
    [notice] Building container 'common/data_processing_dataset_concatenate:dev' with Dockerfile
    [notice] Building container 'label_projection/metrics_accuracy:dev' with Dockerfile
    ...

Viash will build a whole namespace (`ns`) into executables and Nextflow
pipelines into the `target/docker` and `target/nextflow` folders
respectively. By adding the `-q/--query` flag, you can filter which
components to build using a regex. By adding the `--parallel` flag,
these components are built in parallel (otherwise it will take a really
long time). The flag `--setup cachedbuild` will automatically start
building Docker containers for each of these methods.

The command might take a while to run, since it is building a docker
container for each of the components.

**Step 3, run the pipeline with nextflow.** To do so, run the bash
script located at `src/label_projection/workflows/run_nextflow.sh`:

``` bash
src/label_projection/workflows/run/run_test.sh
```

    N E X T F L O W  ~  version 22.04.5
    Launching `src/label_projection/workflows/run/main.nf` [pensive_turing] DSL2 - revision: 16b7b0c332
    executor >  local (28)
    [f6/f89435] process > run_wf:run_methods:true_labels:true_labels_process (pancreas.true_labels)                         [100%] 1 of 1 ✔
    [ed/d674a2] process > run_wf:run_methods:majority_vote:majority_vote_process (pancreas.majority_vote)                   [100%] 1 of 1 ✔
    [15/f0a427] process > run_wf:run_methods:random_labels:random_labels_process (pancreas.random_labels)                   [100%] 1 of 1 ✔
    [02/969d05] process > run_wf:run_methods:knn:knn_process (pancreas.knn)                                                 [100%] 1 of 1 ✔
    [90/5fdf9a] process > run_wf:run_methods:mlp:mlp_process (pancreas.mlp)                                                 [100%] 1 of 1 ✔
    [c7/dee2e5] process > run_wf:run_methods:logistic_regression:logistic_regression_process (pancreas.logistic_regression) [100%] 1 of 1 ✔
    [83/3ba0c9] process > run_wf:run_methods:scanvi:scanvi_process (pancreas.scanvi)                                        [100%] 1 of 1 ✔
    [e3/2c298e] process > run_wf:run_methods:seurat_transferdata:seurat_transferdata_process (pancreas.seurat_transferdata) [100%] 1 of 1 ✔
    [d6/7212ab] process > run_wf:run_methods:xgboost:xgboost_process (pancreas.xgboost)                                     [100%] 1 of 1 ✔
    [b6/7dc1a7] process > run_wf:run_metrics:accuracy:accuracy_process (pancreas.scanvi)                                    [100%] 9 of 9 ✔
    [be/7d4da4] process > run_wf:run_metrics:f1:f1_process (pancreas.scanvi)                                                [100%] 9 of 9 ✔
    [89/dcd77a] process > run_wf:aggregate_results:extract_scores:extract_scores_process (combined)                         [100%] 1 of 1 ✔

## Project structure

High level overview: . ├── bin Helper scripts for building the project
and developing a new component. ├── resources_test Datasets for testing
components. If you don’t have this folder, run **Step 1** above. ├── src
Source files for each component in the pipeline. │ ├── common Common
processing components. │ ├── datasets Components and pipelines for
building the ‘Common datasets’ │ ├── label_projection Source files
related to the ‘Label projection’ task. │ └── … Other tasks. └── target
Executables generated by viash based on the components listed under
`src/`. ├── docker Bash executables which can be used from a terminal.
└── nextflow Nextflow modules which can be used as a standalone pipeline
or as part of a bigger pipeline.

Detailed overview of a task folder (e.g. `src/label_projection`):

    src/label_projection/
    ├── api                    Specs for the components in this task.
    ├── control_methods        Control methods which serve as quality control checks for the benchmark.
    ├── docs                   Task documentation
    ├── methods                Label projection method components.
    ├── metrics                Label projection metric components.
    ├── resources_scripts      The scripts needed to run the benchmark.
    ├── resources_test_scripts The scripts needed to generate the test resources (which are needed for unit testing).
    ├── split_dataset          A component that masks a common dataset for use in the benchmark
    └── workflows              The benchmarking workflow.

Detailed overview of the `src/datasets` folder:

    src/datasets/
    ├── api                    Specs for the data loaders and normalisation methods.
    ├── loaders                Components for ingesting datasets from a source.
    ├── normalization          Normalization method components.
    ├── processors             Other preprocessing components (e.g. HVG and PCA).
    ├── resource_scripts       The scripts needed to generate the common datasets.
    ├── resource_test_scripts  The scripts needed to generate the test resources (which are needed for unit testing).
    └── workflows              The workflow which generates the common datasets.

## Adding a Viash component

[Viash](https://viash.io) allows you to create pipelines in Bash or
Nextflow by wrapping Python, R, or Bash scripts into reusable
components.

You can start creating a new component by [creating a Viash
component](https://viash.io/guide/component/creation/docker.html).

For example, to create a new Python-based method named `foo`, create a
Viash config at `src/label_projection/methods/foo/config.vsh.yaml`:

``` yaml
__merge__: ../../api/comp_method.yaml
functionality:
  name: "foo"
  namespace: "label_projection/methods"
  # A multiline description of your method.
  description: "Todo: fill in"
  info:
    type: method

    # a short label of your method
    label: Foo

    paper_doi: "10.1234/1234.5678.1234567890"

    # if you don't have a Doi, you can specify a name, url and year manually:
    # paper_name: "Nearest neighbor pattern classification"
    # paper_url: "https://doi.org/10.1109/TIT.1967.1053964"
    # paper_year: 1967
    
    code_url: "https://github.com/my_organisation/foo"
  resources:
    - type: python_script
      path: script.py
platforms:
  - type: docker
    image: "python:3.10"
    setup:
      - type: python
        packages:
          - anndata>=0.8
          - scikit-learn
  - type: nextflow
```

And create a script at `src/label_projection/methods/foo/script.py`:

``` python
import anndata as ad
import numpy as np

## VIASH START
# This code-block will automatically be replaced by Viash at runtime.
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'foo'
}
## VIASH END

print("Load data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Create predictions")
input_test.obs["label_pred"] = "foo"

print("Add method name to uns")
input_test.uns["method_id"] = meta["functionality_name"]

print("Write output to file")
input_test.write_h5ad(par["output"], compression="gzip")
```

## Running a component from CLI

You can view the interface of the executable by running the executable
with the `-h` or `--help` parameter.

``` bash
viash run src/label_projection/methods/foo/config.vsh.yaml -- --help
```

    foo dev

    Todo: fill in

    Arguments:
        --input_train
            type: file, file must exist
            example: training.h5ad
            The training data

        --input_test
            type: file, file must exist
            example: test.h5ad
            The test data (without labels)

        --output
            type: file, output, file must exist
            example: prediction.h5ad
            The prediction file

You can **run the component** as follows:

``` bash
viash run src/label_projection/methods/foo/config.vsh.yaml -- \
  --input_train resources_test/label_projection/pancreas/train.h5ad \
  --input_test resources_test/label_projection/pancreas/test.h5ad \
  --output resources_test/label_projection/pancreas/prediction.h5ad
```

    [notice] Checking if Docker image is available at 'ghcr.io/openproblems-bio/label_projection/methods_foo:dev'
    [warning] Could not pull from 'ghcr.io/openproblems-bio/label_projection/methods_foo:dev'. Docker image doesn't exist or is not accessible.
    [notice] Building container 'ghcr.io/openproblems-bio/label_projection/methods_foo:dev' with Dockerfile
    Load data
    Create predictions
    Add method name to uns
    Write output to file

## Building a component

`viash` has several helper functions to help you quickly develop a
component.

With **`viash build`**, you can turn the component into a standalone
executable. This standalone executable you can give to somebody else,
and they will be able to run it, provided that they have Bash and Docker
installed.

``` bash
viash build src/label_projection/methods/foo/config.vsh.yaml \
  -o target/docker/label_projection/methods/foo
```

<div>

> **Note**
>
> The `viash_build` component does a much better job of setting up a
> collection of components.

</div>

You can now view the same interface of the executable by running the
executable with the `-h` parameter.

``` bash
target/docker/label_projection/methods/foo/foo -h
```

    foo dev

    Todo: fill in

    Arguments:
        --input_train
            type: file, file must exist
            example: training.h5ad
            The training data

        --input_test
            type: file, file must exist
            example: test.h5ad
            The test data (without labels)

        --output
            type: file, output, file must exist
            example: prediction.h5ad
            The prediction file

Or **run the component** as follows:

``` bash
target/docker/label_projection/methods/foo/foo \
  --input_train resources_test/label_projection/pancreas/train.h5ad \
  --input_test resources_test/label_projection/pancreas/test.h5ad \
  --output resources_test/label_projection/pancreas/prediction.h5ad
```

    Load data
    Create predictions
    Add method name to uns
    Write output to file

## Unit testing a component

The [method API
specifications](src/label_projection/api/comp_method.yaml) comes with a
generic unit test for free. This means you can unit test your component
using the **`viash test`** command.

``` bash
viash test src/label_projection/methods/foo/config.vsh.yaml
```

    Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo17760509097337858011'
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo17760509097337858011/build_executable/foo ---verbosity 6 ---setup cachedbuild
    [notice] Building container 'ghcr.io/openproblems-bio/label_projection/methods_foo:test_nGvjdE' with Dockerfile
    [info] Running 'docker build -t ghcr.io/openproblems-bio/label_projection/methods_foo:test_nGvjdE /home/rcannood/workspace/viash_temp/viash_test_foo17760509097337858011/build_executable -f /home/rcannood/workspace/viash_temp/viash_test_foo17760509097337858011/build_executable/tmp/dockerbuild-foo-C7VuUU/Dockerfile'
    Sending build context to Docker daemon  39.94kB

    Step 1/7 : FROM python:3.10
     ---> 465483cdaa4e
    Step 2/7 : RUN pip install --upgrade pip &&   pip install --upgrade --no-cache-dir "anndata>=0.8" "scikit-learn"
     ---> Using cache
     ---> 91f658ec0590
    Step 3/7 : LABEL org.opencontainers.image.description="Companion container for running component label_projection/methods foo"
     ---> Using cache
     ---> f1ace85a71b0
    Step 4/7 : LABEL org.opencontainers.image.created="2022-12-17T08:47:34+01:00"
     ---> Running in 299ea3924905
    Removing intermediate container 299ea3924905
     ---> 6fc97da56de8
    Step 5/7 : LABEL org.opencontainers.image.source="https://github.com/openproblems-bio/openproblems-v2"
     ---> Running in bf60068c5fe8
    Removing intermediate container bf60068c5fe8
     ---> 20ff545ec27a
    Step 6/7 : LABEL org.opencontainers.image.revision="8a4877920fc79009dcb1e4bb16674b3b441c75ab"
     ---> Running in c4410d3a7c78
    Removing intermediate container c4410d3a7c78
     ---> 1a57a0d9a7e5
    Step 7/7 : LABEL org.opencontainers.image.version="test_nGvjdE"
     ---> Running in 81d7a66aa40a
    Removing intermediate container 81d7a66aa40a
     ---> 9d84592b1c1e
    Successfully built 9d84592b1c1e
    Successfully tagged ghcr.io/openproblems-bio/label_projection/methods_foo:test_nGvjdE
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo17760509097337858011/test_generic_test/test_executable
    >> Running script as test
    >> Checking whether output file exists
    >> Reading h5ad files
    input_test: AnnData object with n_obs × n_vars = 130 × 443
        obs: 'batch'
        var: 'hvg', 'hvg_score'
        uns: 'dataset_id', 'normalization_id'
        obsm: 'X_pca'
        layers: 'counts', 'normalized'
    output: AnnData object with n_obs × n_vars = 130 × 443
        obs: 'batch', 'label_pred'
        var: 'hvg', 'hvg_score'
        uns: 'dataset_id', 'method_id', 'normalization_id'
        obsm: 'X_pca'
        layers: 'counts', 'normalized'
    >> Checking whether predictions were added
    Checking whether data from input was copied properly to output
    All checks succeeded!
    ====================================================================
    [32mSUCCESS! All 1 out of 1 test scripts succeeded![0m
    Cleaning up temporary directory

Let’s introduce a bug in the script and try running the test again. For
instance:

``` python
import anndata as ad
import numpy as np

## VIASH START
# This code-block will automatically be replaced by Viash at runtime.
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'foo'
}
## VIASH END

print("Load data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Not creating any predictions!!!")
# input_test.obs["label_pred"] = "foo"

print("Not adding method name to uns!!!")
# input_test.uns["method_id"] = meta["functionality_name"]

print("Write output to file")
input_test.write_h5ad(par["output"], compression="gzip")
```

If we now run the test, we should get an error since we didn’t create
all of the required output slots.

``` bash
viash test src/label_projection/methods/foo/config.vsh.yaml
```

    Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128'
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128/build_executable/foo ---verbosity 6 ---setup cachedbuild
    [notice] Building container 'ghcr.io/openproblems-bio/label_projection/methods_foo:test_lnevgh' with Dockerfile
    [info] Running 'docker build -t ghcr.io/openproblems-bio/label_projection/methods_foo:test_lnevgh /home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128/build_executable -f /home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128/build_executable/tmp/dockerbuild-foo-VUcsWQ/Dockerfile'
    Sending build context to Docker daemon  39.94kB

    Step 1/7 : FROM python:3.10
     ---> 465483cdaa4e
    Step 2/7 : RUN pip install --upgrade pip &&   pip install --upgrade --no-cache-dir "anndata>=0.8" "scikit-learn"
     ---> Using cache
     ---> 91f658ec0590
    Step 3/7 : LABEL org.opencontainers.image.description="Companion container for running component label_projection/methods foo"
     ---> Using cache
     ---> f1ace85a71b0
    Step 4/7 : LABEL org.opencontainers.image.created="2022-12-17T08:47:52+01:00"
     ---> Running in ae1e366b6410
    Removing intermediate container ae1e366b6410
     ---> 458c1b49e8b4
    Step 5/7 : LABEL org.opencontainers.image.source="https://github.com/openproblems-bio/openproblems-v2"
     ---> Running in 06a244e7be1e
    Removing intermediate container 06a244e7be1e
     ---> cc48147df9e8
    Step 6/7 : LABEL org.opencontainers.image.revision="8a4877920fc79009dcb1e4bb16674b3b441c75ab"
     ---> Running in 2372d2bddd3d
    Removing intermediate container 2372d2bddd3d
     ---> 7bcad47b5d1b
    Step 7/7 : LABEL org.opencontainers.image.version="test_lnevgh"
     ---> Running in 6499fcfa63af
    Removing intermediate container 6499fcfa63af
     ---> 2213de86e5bc
    Successfully built 2213de86e5bc
    Successfully tagged ghcr.io/openproblems-bio/label_projection/methods_foo:test_lnevgh
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128/test_generic_test/test_executable
    Traceback (most recent call last):
    >> Running script as test
    >> Checking whether output file exists
      File "/viash_automount/home/rcannood/workspace/viash_temp/viash_test_foo4037451094287802128/test_generic_test/tmp//viash-run-foo-BVsf0m.py", line 57, in <module>
    >> Reading h5ad files
        assert "label_pred" in output.obs
    input_test: AnnData object with n_obs × n_vars = 130 × 443
    AssertionError
        obs: 'batch'
        var: 'hvg', 'hvg_score'
        uns: 'dataset_id', 'normalization_id'
        obsm: 'X_pca'
        layers: 'counts', 'normalized'
    output: AnnData object with n_obs × n_vars = 130 × 443
        obs: 'batch'
        var: 'hvg', 'hvg_score'
        uns: 'dataset_id', 'normalization_id'
        obsm: 'X_pca'
        layers: 'counts', 'normalized'
    >> Checking whether predictions were added
    ====================================================================
    [31mERROR! Only 0 out of 1 test scripts succeeded![0m
    Unexpected error occurred! If you think this is a bug, please post
    create an issue at https://github.com/viash-io/viash/issues containing
    a reproducible example and the stack trace below.

    viash - 0.6.6
    Stacktrace:
    java.lang.RuntimeException: Only 0 out of 1 test scripts succeeded!
        at io.viash.ViashTest$.apply(ViashTest.scala:111)
        at io.viash.Main$.internalMain(Main.scala:185)
        at io.viash.Main$.main(Main.scala:77)
        at io.viash.Main.main(Main.scala)

## More information

The [Viash reference docs](https://viash.io/reference/config/) page
provides information on all of the available fields in a Viash config,
and the [Guide](https://viash.io/guide/) will help you get started with
creating components from scratch.

<!-- cleaning up temporary files -->
