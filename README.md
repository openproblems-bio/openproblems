# opsca-viash
Proof Of Concept in adapting Open Problems to use viash.

## Requirements
To use this repository, make sure you have Bash, Java, and Docker installed. If you wish to use Nextflow, install that too.

## Building the pipeline
To build all the components in `src` to `target`, run:

```sh
bin/project_build
```

Note that this will also build docker containers for each of the components. The first time, this might take a while.

## Running a simple pipeline with Bash

Inspect the contents of the sample Bash pipeline in `run_bash_pipeline.sh`, then run it.

```sh
./run_bash_pipeline.sh
```

Some components might take a while to run. In the end, the `out_bash/` folder will have been populated with a lot of h5ad files, and a `scores.tsv` file.


## Running a simple pipeline with Nextflow

Inspect the contents of the sample bash pipeline in `run_nxf_pipeline.sh`, then run it.

```
$ ./run_nxf_pipeline.sh
N E X T F L O W  ~  version 20.04.1-edge
Launching `main.nf` [suspicious_lavoisier] - revision: f6dca0e8d5
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
executor >  local (3)
[16/e727b4] process > citeseq_cbmc:citeseq_cbmc_process (dyngen) [100%] 1 of 1 ✔
[85/88774c] process > mnn:mnn_process (dyngen)                   [100%] 1 of 1 ✔
[1e/9593dd] process > knn_auc:knn_auc_process (dyngen)           [100%] 1 of 1 ✔
WARN: Access to undefined parameter `debug` -- Initialise it to a default value eg. `params.debug = some_value`
Completed at: 06-Apr-2021 22:52:43
Duration    : 2m 46s
CPU hours   : (a few seconds)
Succeeded   : 3
```

Again, some components might take a while to run. 


## Component development with viash

You can run, build, or test a component individually with `bin/viash`.

A tutorial on how to create components with viash can be found at
[github.com/data-intuitive/viash\_tutorial\_1](https://github.com/data-intuitive/viash_tutorial_1).

More documentation is available at
[data-intuitive.com/viash\_docs](https://www.data-intuitive.com/viash_docs) (WIP).

### Creating a new component

<!-- todo: expand documentation -->

Create a viash config file and write a script.

### View help of a component

``` bash
viash run src/modality_alignment/methods/mnn/config.vsh.yaml -- -h
```

Or if you’ve already built the component (see below)

``` bash
target/docker/modality_alignment/methods/mnn/mnn -h
```

### Build a component

``` bash
viash build src/modality_alignment/methods/mnn/config.vsh.yaml -p docker -o target/docker/modality_alignment/methods/mnn --setup
```

### Test a component

``` bash
viash test src/modality_alignment/methods/mnn/config.vsh.yaml
```

### Run a component

``` bash
viash run src/modality_alignment/methods/mnn/config.vsh.yaml -- [... arguments for the component ...]
```

Or if you’ve already built the component

``` bash
target/docker/modality_alignment/methods/mnn/mnn [... arguments for the component ...]
```
