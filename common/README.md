# common_resources
This repo contains common resources that can be used for OpenProblems v2 tasks.

## Usage

> [!NOTE]
> The following instructions are not required when using the [task_template](https://github.com/openproblems-bio/task_template) repository to create your task repository.

To use the resources in this repository, you will need to add this as a submodule to the task repository. 

You can do this by running the following command:

```bash
git submodule add git@github.com:openproblems-bio/common_resources.git common
```

## Update

To update the repository with the latest changes from in the submodule, you can run the following command:

```bash
git submodule update --remote
```

## Initialize

When cloning a repository with a submodule and there are no files visible, you will need to initialize by running the following command:

```bash
git submodule update --init --recursive
```

## Resources

The above information is also available on [working-with-submodules](https://github.blog/2016-02-01-working-with-submodules/).
