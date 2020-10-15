<!--- TODO: update --->

# Template Dataset

Here's a brief task description, maybe link to some seminal papers.

## API

Datasets should contain the following attributes:

* `adata.obs["template_variable"]`

Methods should assign output to `adata.obs['template_output']`.

Metrics should compare `adata['template_output']` to `adata.obs['template_variable']`.
