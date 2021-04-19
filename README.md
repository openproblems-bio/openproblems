# opsca-viash
Proof Of Concept in adapting Open Problems to use viash.

## Requirements
To use this repository, make sure you have Bash, Java, and Docker installed. If you wish to use Nextflow, install that too.

## Quick start

Running the modality alignment pipeline requires two simple steps.

First, by running the command below, viash will build all the components in the `src/` folder as executables in the `target/` folder.

```bash
bin/project_build
```

Next, to run the pipeline with nextflow, run the bash script located at `src/modality_alignment/workflows/run_nextflow.sh`:

```
$ src/modality_alignment/workflows/run_nextflow.sh 
[86/3d1927] process > scprep_csv:scprep_csv_process (CBMC_8K_13AB_10x) [100%] 1 of 1 ✔
[00/b3528e] process > mnn:mnn_process (CBMC_8K_13AB_10x)               [100%] 1 of 1 ✔
[75/bfdbb1] process > scot:scot_process (CBMC_8K_13AB_10x)             [100%] 1 of 1 ✔
[bc/da3dce] process > knn_auc:knn_auc_process (CBMC_8K_13AB_10x.scot)  [100%] 2 of 2 ✔
[b9/e083fe] process > extract_scores:extract_scores_process (combined) [100%] 1 of 1 ✔
Completed at: 19-Apr-2021 13:18:26
Duration    : 3m 30s
CPU hours   : 0.1
Succeeded   : 6
```

## Project structure

```
bin/                     Helper scripts for building the project and developing a new component.
src/                     Source files for each component in the pipeline.
  modality_alignment/    Source files related to the 'Modality alignment' task.
    datasets/            Dataset downloader components.
    methods/             Modality alignment method components.
    metrics/             Modality alignment metric components.
    resources/           Helper files.
    workflow/            The pipeline workflow for this task.
  utils/                 Helper files.
target/                  Executables generated by viash based on the components listed under `src/`.
  docker/                Bash executables which can be used from a terminal.
  nextflow/              Nextflow modules which can be used in a Nextflow pipeline.
work/                    A working directory used by Nextflow.
output/                  Output generated by the pipeline.
```


## Adding a component with viash

[`viash`](https://github.com/data-intuitive/viash) allows you to create pipelines
in Bash or Nextflow by wrapping Python, R, or Bash scripts into reusable components.

### Create a skeleton component
You can start creating a new component using `bin/skeleton`. For example to create
a new Python-based viahs component in the `src/modality_alignment/methods/foo` folder, run:
You can start creating a new component by using the `bin/skeleton` command:

```bash
bin/skeleton --name foo --namespace "modality_alignment/methods" --language python

# or:
bin/skeleton -n foo -ns "modality_alignment/methods" -l python
```

This should create a few files in this folder:

```
script.py                A python script for you to edit.
config.vsh.yaml          Metadata for the script containing info on the input/output arguments of the component.
test.py                  A python script with which you can start unit testing your component.
```

The [Getting started](http://www.data-intuitive.com/viash_docs/) page on the viash documentation site
provides some information on how a basic viash component works, or on the specifications of the `config.vsh.yaml` [config file](http://www.data-intuitive.com/viash_docs/config/).

### Building a component

`viash` has several helper functions to help you quickly develop a component.

With **`viash build`**, you can turn the component into a standalone executable. 
This standalone executable you can give to somebody else, and they will be able to 
run it, provided that they have Bash and Docker installed.

```
viash build src/modality_alignment/methods/foo/config.vsh.yaml \
  -o target/docker/modality_alignment/methods/foo \
  --setup
```

Note that the `bin/project_build` component does a much better job of setting up 
a collection of components. You can filter which components will be built by 
providing a regex to the `-q` parameter, e.g. `bin/project_build -q 'utils|modality_alignment'`.

### Running a component from CLI

You can view the interface of the executable by running the executable with the `-h` parameter.

```
$ target/docker/modality_alignment/methods/foo/foo -h
Replace this with a (multiline) description of your component.

Options:
    -i file, --input=file
        type: file, required parameter
        Describe the input file.

    -o file, --output=file
        type: file, required parameter
        Describe the output file.

    --option=string
        type: string, default: default-
        Describe an optional parameter.
```

You can **run the component** as follows:

```
$ target/docker/modality_alignment/methods/foo/foo -i LICENSE -o output.txt
This is a skeleton component
The arguments are:
 - input:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/LICENSE
 - output:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/output.txt
 - option:  default-
```

Alternatively, you can run the component straight from the viash config by using the **`viash run`** command:
```
viash build src/modality_alignment/methods/foo/config.vsh.yaml -- -i LICENSE -o output.txt
```

### Unit testing a component
Provided that you wrote a script that allows you to test the functionality of a component, 
you can run the tests by using the **`viash test`** command.

```
$ viash test src/modality_alignment/methods/foo/config.vsh.yaml
Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo8028146580425979678'
====================================================================
+/home/rcannood/workspace/viash_temp/viash_test_foo8028146580425979678/build_executable/foo ---setup
> docker build -t modality_alignment/methods_foo:0.0.1 /home/rcannood/workspace/viash_temp/viashsetupdocker-foo-RNMkfg
====================================================================
+/home/rcannood/workspace/viash_temp/viash_test_foo8028146580425979678/test_test.py/test.py
.
----------------------------------------------------------------------
Ran 1 test in 0.016s

OK
====================================================================
SUCCESS! All 1 out of 1 test scripts succeeded!
Cleaning up temporary directory
```

To run all the unit tests of all the components in the repository, use `bin/project_test`. 

## Frequently asked questions

### Running a component causes error 'Unable to find image'

Depending on how an executable was created, a Docker container might not have been created. 

To solve this issue, run the executable with a `---setup` flag attached. This will 
automatically build the Docker container for you. 

```
$ target/docker/modality_alignment/methods/foo/foo ---setup
> docker build -t modality_alignment/methods_foo:0.0.1 /home/rcannood/workspace/viash_temp/viashsetupdocker-foo-KeBjFs
```

Or when working with `viash run`:

```
$ viash run src/modality_alignment/methods/mnn/config.vsh.yaml -- ---setup
```
