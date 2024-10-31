# Instructions

This is a guide on what to do after you have created a new task repository from the template. More in depth information about how to create a new task can be found in the [OpenProblems Documentation](https://openproblems.bio/documentation/create_task/).

## Requirements

A list of required software to start developing a new task can be found in the [OpenProblems Requirements](https://openproblems.bio/documentation/create_task/requirements).

## First things first

### `_viash.yaml`

1. Update the `name` field to the name of the task in snake_case the name should start with `task_`.
2. Update the `description` field to a short description of the task.
3. Add a keyword to the `keywords` field that describes the task.
4. Update the `task_template` in the links fields to the name of the task from step 1.
5. Update the info fields to the text from the task issue.
6. Replace the task_template to the name of the task in `info.name`.
7. Update the authors of the task.

After performing the steps above you can remove the comments in the `_viash.yaml` file.

### `common` submodule

Check the instructions [here](README.md) for more information.

## Resources

The OpenProblems team has provided some test resources that can be used to test the task. These resources are stored in the `resources_test` folder. The `scripts/download_resources.sh` script can be used to download these resources.
If these resources are not sufficient, you can add more resources to the `resources_test` folder. The `scripts/download_resources.sh` script can be updated to download these resources. When using new test resources let the OP team know so they can be added to the s3 bucket.

```bash	
scripts/download_resources.sh
```

## Next steps

### API files ([docs](https://openproblems.bio/documentation/create_task/design_api))

Update the API files in the `src/api` folder. These files define the input and output of the methods and metrics. 

### Components ([docs](https://openproblems.bio/documentation/create_task/create_components))

To create a component, you can run the respective script in the `script` directory. Before running the script make sure to update the variables `task_name`, `component_name` and `component_lang` and save the file. For additional components you will only need to update the `component_name` and `component_lang` variables.

```bash
scripts/add_a_control_method.sh
```

```bash
scripts/add_a_method.sh
```

```bash
scripts/add_a_metric.sh
```

For each type of component there already is a first component created that you can modify.

1. Update the `.info` fields in the  `config.vsh.yaml`.
2. Add any component specific arguments to the `config.vsh.yaml` file.
3. Add any additional resources that are required for the component.
4. Update the docker engine image setup if additional packages are required.
5. If you know the required memory and or CPU you can adjust the nextflow `.directive.labels` field. In addition if your component requires a GPU you can add the `gpu` label to the field.
6. Update the `script.py` or `script.R` file with the code for the component.

> [!NOTE]
> You can remove the comments in the `config.vsh.yaml` file after you have updated the file.

### Testing Components ([docs](https://openproblems.bio/documentation/create_component/run_tests))

You can test the component by running the following command:

```bash
viash test /path/to/config.vsh.yaml
```

Y0u can also test all components by running the following command:

```bash
scripts/test_all_components.sh
```

It is possible to customise the command in the above script by adding a `-q` argument to only perform the test on for example methods e.g. ` -q methods`.


## Dataset processor ([docs](https://openproblems.bio/documentation/create_task/dataset_processor))

The dataset processor is a script that removes all unnecessary info from the dataset for your task. This info is defined in the `api/file_common_dataset.yaml` file. From this filtered dataset several files are created that are used by the methods and metrics. Safeguarding data leaks and laking sure the structure of the data cannot be altered for a method or a metric.

To create the data processor there is no template available. You can follow the guideline in the documentation. Store the processor in the `src/process_dataset` folder.

Be sure to update the `file_common_dataset.yaml` file with the correct information required for the methods/metrics.

> [!IMPORTANT]
> When using your own datasets please advise the OpenProblems team on how to add these datasets to the s3 bucket.
> As the dataset processor should make use of the `common` datasets folder in the `resources` or `resources_test` directory.

To create the resources and test_resources for the task we will create a nextflow workflow that will process the datasets. This workflow will be created together with the OpenProblems team.

## README

To create the task `README` file preform following command:

```bash
scripts/create_readme.sh
```

## Benchmarking ([docs](https://openproblems.bio/documentation/create_task/create_workflow))

When you are finished with creating your components and dataset processor you can create a workflow to benchmark the components. This workflow will be created together with the OpenProblems team.
