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

**Step 0, fetch viash and nextflow:** run the `bin/init` executable.

``` bash
bin/init
```

    > Using tag develop
    > Cleanup
    > Downloading Viash source code @develop
    > Building Viash from source
    > Building Viash helper scripts from source
    > Done, happy viash-ing!

**Step 1, download test resources:** by running the following command.

``` bash
bin/viash run src/common/sync_test_resources/config.vsh.yaml
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
bin/viash_build -q 'label_projection|common'
```

    In development mode with 'dev'.
    Exporting split_dataset (label_projection) =docker=> target/docker/label_projection/split_dataset
    Exporting accuracy (label_projection/metrics) =docker=> target/docker/label_projection/metrics/accuracy
    Exporting random_labels (label_projection/control_methods) =docker=> target/docker/label_projection/control_methods/random_labels
    [notice] Building container 'label_projection/control_methods_random_labels:dev' with Dockerfile
    [notice] Building container 'common/data_processing_dataset_concatenate:dev' with Dockerfile
    [notice] Building container 'label_projection/metrics_accuracy:dev' with Dockerfile
    ...

These standalone executables you can give to somebody else, and they
will be able to run it, provided that they have Bash and Docker
installed. The command might take a while to run, since it is building a
docker container for each of the components.

**Step 3, run the pipeline with nextflow.** To do so, run the bash
script located at `src/label_projection/workflows/run_nextflow.sh`:

``` bash
src/label_projection/workflows/run/run_test.sh
```

    N E X T F L O W  ~  version 22.04.5
    Launching `src/label_projection/workflows/run/main.nf` [small_becquerel] DSL2 - revision: ece87259df
    executor >  local (19)
    [39/e1bb01] process > run_wf:true_labels:true_labels_process (1)                 [100%] 1 of 1 ✔
    [3b/d41f8a] process > run_wf:random_labels:random_labels_process (1)             [100%] 1 of 1 ✔
    [c2/0398dd] process > run_wf:majority_vote:majority_vote_process (1)             [100%] 1 of 1 ✔
    [fd/92edc7] process > run_wf:knn_classifier:knn_classifier_process (1)           [100%] 1 of 1 ✔
    [f7/7cdb34] process > run_wf:logistic_regression:logistic_regression_process (1) [100%] 1 of 1 ✔
    [4f/6a67e4] process > run_wf:mlp:mlp_process (1)                                 [100%] 1 of 1 ✔
    [a5/ae6341] process > run_wf:accuracy:accuracy_process (6)                       [100%] 6 of 6 ✔
    [72/5076e8] process > run_wf:f1:f1_process (6)                                   [100%] 6 of 6 ✔
    [cf/eccd48] process > run_wf:extract_scores:extract_scores_process               [100%] 1 of 1 ✔

## Project structure

    .
    ├── bin                    Helper scripts for building the project and developing a new component.
    ├── resources_test         Datasets for testing components. If you don't have this folder, run **Step 1** above.
    ├── src                    Source files for each component in the pipeline.
    │   ├── common             Common processing components.
    │   ├── datasets           Components for ingesting datasets from a source.
    │   ├── label_projection   Source files related to the 'Label projection' task.
    │   └── ...                Other tasks.
    └── target                 Executables generated by viash based on the components listed under `src/`.
        ├── docker             Bash executables which can be used from a terminal.
        └── nextflow           Nextflow modules which can be used as a standalone pipeline or as part of a bigger pipeline.


    bin/                     Helper scripts for building the project and developing a new component.
    resources_test/          Datasets for testing components.
    src/                     Source files for each component in the pipeline.
      common/                Common processing components.
      datasets/              Components related to ingesting datasets into OpenProblems v2.
        api/                 Specs for the data loaders and normalisation methods.
        loaders/             Components for ingesting datasets from a source.
        normalization/       Common normalization methods.
      label_projection/      Source files related to the 'Label projection' task.
        datasets/            Dataset downloader components.
        methods/             Modality alignment method components.
        metrics/             Modality alignment metric components.
        utils/               Utils functions.
        workflow/            The pipeline workflow for this task.
    target/                  Executables generated by viash based on the components listed under `src/`.
      docker/                Bash executables which can be used from a terminal.
      nextflow/              Nextflow modules which can be used in a Nextflow pipeline.
    work/                    A working directory used by Nextflow.
    output/                  Output generated by the pipeline.

The `src/datasets` folder

src/datasets/ ├── api Specs for the data loaders and normalisation
methods. ├── loaders Components for ingesting datasets from a source.
├── normalization Common normalization methods. ├──
resource_test_scripts Scripts for generating the objects in the
`resources_test` folder. └── workflows A set of Nextflow workflows which
tie together various components.

The `src/label_projection` folder

src/label_projection/ ├── api Specs for the split_dataset, methods and
metrics in this task. ├── control_methods Positive and negative control
methods for quality control. ├── methods Method components. ├── metrics
Metric components. ├── [README.md](src/label_projection/) More
information on how this task works. ├── resources_test_scripts Scripts
for generating the objects in the `resources_test` folder. ├──
split_dataset A component for splitting a common dataset into a `train`,
`test` and `solution` object. └── workflows A set of Nextflow workflows
which tie together various components.

## Adding a Viash component

[Viash](https://viash.io) allows you to create pipelines in Bash or
Nextflow by wrapping Python, R, or Bash scripts into reusable
components.

You can start creating a new component by [creating a Viash
component](https://viash.io/guide/component/creation/docker.html).

For example, to create a new Python-based method named `foo`, create a
Viash config at `src/label_projection/methods/foo/config.vsh.yaml`:

``` yaml
__inherits__: ../../api/comp_method.yaml
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
bin/viash run src/label_projection/methods/foo/config.vsh.yaml -- --help
```

    Warning: Config inheritance (__inherits__) is an experimental feature. Changes to the API are expected.
    foo

    Todo: fill in

    Arguments:
        --input_train
            type: file
            example: training.h5ad
            The training data

        --input_test
            type: file
            example: test.h5ad
            The test data (without labels)

        --output
            type: file, output
            example: prediction.h5ad
            The prediction file

You can **run the component** as follows:

``` bash
bin/viash run src/label_projection/methods/foo/config.vsh.yaml -- \
  --input_train resources_test/label_projection/pancreas/train.h5ad \
  --input_test resources_test/label_projection/pancreas/test.h5ad \
  --output resources_test/label_projection/pancreas/prediction.h5ad
```

    Warning: Config inheritance (__inherits__) is an experimental feature. Changes to the API are expected.
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
bin/viash build src/label_projection/methods/foo/config.vsh.yaml \
  -o target/docker/label_projection/methods/foo
```

    Warning: Config inheritance (__inherits__) is an experimental feature. Changes to the API are expected.

<div>

> **Note**
>
> The `bin/viash_build` component does a much better job of setting up a
> collection of components.

</div>

You can now view the same interface of the executable by running the
executable with the `-h` parameter.

``` bash
target/docker/label_projection/methods/foo/foo -h
```

    foo

    Todo: fill in

    Arguments:
        --input_train
            type: file
            example: training.h5ad
            The training data

        --input_test
            type: file
            example: test.h5ad
            The test data (without labels)

        --output
            type: file, output
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
bin/viash test src/label_projection/methods/foo/config.vsh.yaml
```

    Warning: Config inheritance (__inherits__) is an experimental feature. Changes to the API are expected.
    Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo7865291233056269818'
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo7865291233056269818/build_executable/foo ---verbosity 6 ---setup cachedbuild
    [notice] Building container 'label_projection/methods_foo:test_rIrBSI' with Dockerfile
    [info] Running 'docker build -t label_projection/methods_foo:test_rIrBSI /home/rcannood/workspace/viash_temp/viash_test_foo7865291233056269818/build_executable -f /home/rcannood/workspace/viash_temp/viash_test_foo7865291233056269818/build_executable/tmp/dockerbuild-foo-y7Hdos/Dockerfile'
    Sending build context to Docker daemon  37.89kB

    Step 1/7 : FROM python:3.10
     ---> ecbdd6bafdb5
    Step 2/7 : RUN pip install --upgrade pip &&   pip install --upgrade --no-cache-dir "anndata>=0.8" "scikit-learn"
     ---> Using cache
     ---> f1fbd09c8ccd
    Step 3/7 : LABEL org.opencontainers.image.description="Companion container for running component label_projection/methods foo"
     ---> Using cache
     ---> 063049300b14
    Step 4/7 : LABEL org.opencontainers.image.created="2022-11-17T07:37:39+01:00"
     ---> Running in 5d3bda2ec79c
    Removing intermediate container 5d3bda2ec79c
     ---> e77cbb1b5502
    Step 5/7 : LABEL org.opencontainers.image.source="https://github.com/openproblems-bio/openproblems-v2.git"
     ---> Running in fed5d3371cea
    Removing intermediate container fed5d3371cea
     ---> 9db7959fd7af
    Step 6/7 : LABEL org.opencontainers.image.revision="8f9371ddfa5f5c20df01612342040a2003274da3"
     ---> Running in 9a2f654aedb9
    Removing intermediate container 9a2f654aedb9
     ---> 912628903c90
    Step 7/7 : LABEL org.opencontainers.image.version="test_rIrBSI"
     ---> Running in 82b6f860d949
    Removing intermediate container 82b6f860d949
     ---> e0720a8407a9
    Successfully built e0720a8407a9
    Successfully tagged label_projection/methods_foo:test_rIrBSI
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo7865291233056269818/test_generic_test/test_executable
    >> Running script as test
    >> Checking whether output file exists
    >> Reading h5ad files
    input_test: AnnData object with n_obs × n_vars = 307 × 443
        obs: 'batch'
        uns: 'dataset_id'
        layers: 'counts', 'log_cpm', 'log_scran_pooling'
    output: AnnData object with n_obs × n_vars = 307 × 443
        obs: 'batch', 'label_pred'
        uns: 'dataset_id', 'method_id'
        layers: 'counts', 'log_cpm', 'log_scran_pooling'
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
bin/viash test src/label_projection/methods/foo/config.vsh.yaml
```

    Warning: Config inheritance (__inherits__) is an experimental feature. Changes to the API are expected.
    Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789'
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789/build_executable/foo ---verbosity 6 ---setup cachedbuild
    [notice] Building container 'label_projection/methods_foo:test_5q5NGA' with Dockerfile
    [info] Running 'docker build -t label_projection/methods_foo:test_5q5NGA /home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789/build_executable -f /home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789/build_executable/tmp/dockerbuild-foo-GRvLt9/Dockerfile'
    Sending build context to Docker daemon  37.89kB

    Step 1/7 : FROM python:3.10
     ---> ecbdd6bafdb5
    Step 2/7 : RUN pip install --upgrade pip &&   pip install --upgrade --no-cache-dir "anndata>=0.8" "scikit-learn"
     ---> Using cache
     ---> f1fbd09c8ccd
    Step 3/7 : LABEL org.opencontainers.image.description="Companion container for running component label_projection/methods foo"
     ---> Using cache
     ---> 063049300b14
    Step 4/7 : LABEL org.opencontainers.image.created="2022-11-17T07:38:03+01:00"
     ---> Running in 2211c2f3d253
    Removing intermediate container 2211c2f3d253
     ---> 4a7607ecb7b4
    Step 5/7 : LABEL org.opencontainers.image.source="https://github.com/openproblems-bio/openproblems-v2.git"
     ---> Running in 21c85f7d64bc
    Removing intermediate container 21c85f7d64bc
     ---> ad15f8b03066
    Step 6/7 : LABEL org.opencontainers.image.revision="8f9371ddfa5f5c20df01612342040a2003274da3"
     ---> Running in 94c15d5c9736
    Removing intermediate container 94c15d5c9736
     ---> 423d9a6d04e2
    Step 7/7 : LABEL org.opencontainers.image.version="test_5q5NGA"
     ---> Running in d31d99d449ef
    Removing intermediate container d31d99d449ef
     ---> 7f9635195fe8
    Successfully built 7f9635195fe8
    Successfully tagged label_projection/methods_foo:test_5q5NGA
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789/test_generic_test/test_executable
    Traceback (most recent call last):
    >> Running script as test
      File "/viash_automount/home/rcannood/workspace/viash_temp/viash_test_foo15779522933199950789/test_generic_test/tmp//viash-run-foo-j6Jfba", line 56, in <module>
    >> Checking whether output file exists
        assert "label_pred" in output.obs
    >> Reading h5ad files
    AssertionError
    input_test: AnnData object with n_obs × n_vars = 307 × 443
        obs: 'batch'
        uns: 'dataset_id'
        layers: 'counts', 'log_cpm', 'log_scran_pooling'
    output: AnnData object with n_obs × n_vars = 307 × 443
        obs: 'batch'
        uns: 'dataset_id'
        layers: 'counts', 'log_cpm', 'log_scran_pooling'
    >> Checking whether predictions were added
    ====================================================================
    [31mERROR! Only 0 out of 1 test scripts succeeded![0m
    Unexpected error occurred! If you think this is a bug, please post
    create an issue at https://github.com/viash-io/viash/issues containing
    a reproducible example and the stack trace below.

    viash - 0.6.3
    Stacktrace:
    java.lang.RuntimeException: Only 0 out of 1 test scripts succeeded!
        at io.viash.ViashTest$.apply(ViashTest.scala:110)
        at io.viash.Main$.internalMain(Main.scala:99)
        at io.viash.Main$.main(Main.scala:39)
        at io.viash.Main.main(Main.scala)

## More information

The [Viash reference docs](https://viash.io/reference/config/) page
provides information on all of the available fields in a Viash config,
and the [Guide](https://viash.io/guide/) will help you get started with
creating components from scratch.

<!-- cleaning up temporary files -->
