+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 30  # Order that this section will appear.

title = "Label projection"
subtitle = "Using cell labels from a reference dataset to annotate an unseen dataset"

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"


  # Text color (true=light or false=dark).
  text_color_light = false

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS.
 css_style = ""

 # CSS class.
 css_class = ""

+++
|    | metric   |      value | method                      | task             |   memory_mb |   memory_leaked_mb |   runtime_s | dataset          |
|---:|:---------|-----------:|:----------------------------|:-----------------|------------:|-------------------:|------------:|:-----------------|
|  0 | accuracy | 0.554878   | logistic_regression_scran   | Label Projection |   633.137   |        627.602     |  14.9918    | pancreas_batch   |
|  0 | accuracy | 0.567073   | mlp_scran                   | Label Projection |    30.9023  |         24.9648    |   0.586668  | pancreas_batch   |
|  0 | accuracy | 0.579268   | mlp_log_cpm                 | Label Projection |     6.24609 |          0.285156  |   0.579627  | pancreas_batch   |
|  0 | accuracy | 0.621951   | logistic_regression_log_cpm | Label Projection |     7.26562 |          1.48828   |   0.0520863 | pancreas_batch   |
|  1 | f1       | 0.560702   | logistic_regression_scran   | Label Projection |   633.137   |        627.602     |  14.9918    | pancreas_batch   |
|  1 | f1       | 0.572338   | mlp_scran                   | Label Projection |    30.9023  |         24.9648    |   0.586668  | pancreas_batch   |
|  1 | f1       | 0.583483   | mlp_log_cpm                 | Label Projection |     6.24609 |          0.285156  |   0.579627  | pancreas_batch   |
|  1 | f1       | 0.624757   | logistic_regression_log_cpm | Label Projection |     7.26562 |          1.48828   |   0.0520863 | pancreas_batch   |
|  0 | accuracy | 0.693333   | logistic_regression_scran   | Label Projection |     6.13672 |          0.117188  |   0.0792254 | pancreas_random  |
|  0 | accuracy | 0.706667   | mlp_scran                   | Label Projection |    26.4219  |         20.4375    |   0.626644  | pancreas_random  |
|  0 | accuracy | 0.746667   | logistic_regression_log_cpm | Label Projection |     5.98047 |          0.0546875 |   0.0506659 | pancreas_random  |
|  0 | accuracy | 0.8        | mlp_log_cpm                 | Label Projection |     5.96875 |          0         |   0.620226  | pancreas_random  |
|  1 | f1       | 0.697304   | logistic_regression_scran   | Label Projection |     6.13672 |          0.117188  |   0.0792254 | pancreas_random  |
|  1 | f1       | 0.710184   | mlp_scran                   | Label Projection |    26.4219  |         20.4375    |   0.626644  | pancreas_random  |
|  1 | f1       | 0.749658   | logistic_regression_log_cpm | Label Projection |     5.98047 |          0.0546875 |   0.0506659 | pancreas_random  |
|  1 | f1       | 0.801547   | mlp_log_cpm                 | Label Projection |     5.96875 |          0         |   0.620226  | pancreas_random  |
|  0 | accuracy | 0          | logistic_regression_log_cpm | Label Projection |     5.83984 |          0.0976562 |   0.0701691 | zebrafish_labels |
|  0 | accuracy | 0          | logistic_regression_scran   | Label Projection |     5.79688 |          0.0703125 |   0.288389  | zebrafish_labels |
|  0 | accuracy | 0          | mlp_scran                   | Label Projection |     5.74219 |          0         |   0.211888  | zebrafish_labels |
|  0 | accuracy | 0.0163934  | mlp_log_cpm                 | Label Projection |     5.76562 |          0         |   0.244471  | zebrafish_labels |
|  1 | f1       | 0          | logistic_regression_log_cpm | Label Projection |     5.83984 |          0.0976562 |   0.0701691 | zebrafish_labels |
|  1 | f1       | 0          | logistic_regression_scran   | Label Projection |     5.79688 |          0.0703125 |   0.288389  | zebrafish_labels |
|  1 | f1       | 0          | mlp_scran                   | Label Projection |     5.74219 |          0         |   0.211888  | zebrafish_labels |
|  1 | f1       | 0.00252207 | mlp_log_cpm                 | Label Projection |     5.76562 |          0         |   0.244471  | zebrafish_labels |
|  0 | accuracy | 0.115385   | mlp_log_cpm                 | Label Projection |     5.63672 |          0         |   0.208209  | zebrafish_random |
|  0 | accuracy | 0.115385   | mlp_scran                   | Label Projection |     5.62109 |          0         |   0.22847   | zebrafish_random |
|  0 | accuracy | 0.134615   | logistic_regression_log_cpm | Label Projection |     5.71484 |          0         |   0.0692928 | zebrafish_random |
|  0 | accuracy | 0.192308   | logistic_regression_scran   | Label Projection |     5.66406 |          0         |   0.0699402 | zebrafish_random |
|  1 | f1       | 0.0779093  | mlp_scran                   | Label Projection |     5.62109 |          0         |   0.22847   | zebrafish_random |
|  1 | f1       | 0.107166   | mlp_log_cpm                 | Label Projection |     5.63672 |          0         |   0.208209  | zebrafish_random |
|  1 | f1       | 0.136126   | logistic_regression_log_cpm | Label Projection |     5.71484 |          0         |   0.0692928 | zebrafish_random |
|  1 | f1       | 0.154773   | logistic_regression_scran   | Label Projection |     5.66406 |          0         |   0.0699402 | zebrafish_random |