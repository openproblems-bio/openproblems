OpenProblems v2
================

Formalizing and benchmarking open problems in single-cell genomics.

[**Visit the Open Problems Website.**](https://openproblems.bio/)

To get started with developing a new method or metric, please see the
[Contribution guidelines](CONTRIBUTING.md).

## Benefits of using Viash

### The pipeline is **language-agnostic**

This means that each component can be written in whatever scripting
language the user desires. By default, Viash supports wrapping the
following scripting languages: Bash, Python, R, JavaScript, and Scala.
If Viash doesn’t support your preferred scripting language, you can
still write a Bash script which calls your desired programming language.

### One Docker container per component

Viash builds one Docker container per component. While this results in
some initial computational overhead, this makes it a lot easier to add a
new component to the pipeline with dependencies which might conflict
with those of other components.

### Reproducible components

A component built by viash is meant to be reproducible. All executables
and Nextflow modules in the `target/` folder in one of the
[releases](https://github.com/openproblems-bio/openproblems-v2/releases)
is fully reproducible, since all containers are published on the [GitHub
Container
Registry](https://github.com/orgs/openproblems-bio/packages?repo_name=openproblems-v2).
