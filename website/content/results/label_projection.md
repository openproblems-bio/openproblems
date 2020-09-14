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
|    | metric   |    value | method                      | task             |   memory_mb |   memory_leaked_mb |   runtime_s | dataset          |
|---:|:---------|---------:|:----------------------------|:-----------------|------------:|-------------------:|------------:|:-----------------|
|  0 | accuracy | 0.957263 | logistic_regression_log_cpm | Label Projection |     6478.39 |          1236.16   |     16.3349 | pancreas_batch   |
|  0 | accuracy | 0.964295 | mlp_log_cpm                 | Label Projection |     9837.72 |          1193.29   |     19.7333 | pancreas_batch   |
|  0 | accuracy | 0.966189 | logistic_regression_scran   | Label Projection |     9816.19 |          2915.69   |    277.899  | pancreas_batch   |
|  0 | accuracy | 0.969164 | mlp_scran                   | Label Projection |    10844.3  |          1489.36   |    276.515  | pancreas_batch   |
|  1 | f1       | 0.957265 | logistic_regression_log_cpm | Label Projection |     6478.39 |          1236.16   |     16.3349 | pancreas_batch   |
|  1 | f1       | 0.964247 | mlp_log_cpm                 | Label Projection |     9837.72 |          1193.29   |     19.7333 | pancreas_batch   |
|  1 | f1       | 0.96668  | logistic_regression_scran   | Label Projection |     9816.19 |          2915.69   |    277.899  | pancreas_batch   |
|  1 | f1       | 0.969669 | mlp_scran                   | Label Projection |    10844.3  |          1489.36   |    276.515  | pancreas_batch   |
|  0 | accuracy | 0.986343 | logistic_regression_log_cpm | Label Projection |    12000.1  |          1193.79   |     17.3775 | pancreas_random  |
|  0 | accuracy | 0.986343 | logistic_regression_scran   | Label Projection |    11605.1  |          -361.043  |    398.035  | pancreas_random  |
|  0 | accuracy | 0.987557 | mlp_log_cpm                 | Label Projection |    11817.6  |          1193.17   |     37.9832 | pancreas_random  |
|  0 | accuracy | 0.988467 | mlp_scran                   | Label Projection |    11554    |           144.898  |    468.924  | pancreas_random  |
|  1 | f1       | 0.985882 | logistic_regression_log_cpm | Label Projection |    12000.1  |          1193.79   |     17.3775 | pancreas_random  |
|  1 | f1       | 0.986007 | logistic_regression_scran   | Label Projection |    11605.1  |          -361.043  |    398.035  | pancreas_random  |
|  1 | f1       | 0.987224 | mlp_log_cpm                 | Label Projection |    11817.6  |          1193.17   |     37.9832 | pancreas_random  |
|  1 | f1       | 0.988402 | mlp_scran                   | Label Projection |    11554    |           144.898  |    468.924  | pancreas_random  |
|  0 | accuracy | 0.211626 | mlp_scran                   | Label Projection |    13842.6  |           988.41   |   1126.67   | zebrafish_labels |
|  0 | accuracy | 0.220384 | mlp_log_cpm                 | Label Projection |     6325.73 |           376.273  |     89.4715 | zebrafish_labels |
|  0 | accuracy | 0.254249 | logistic_regression_log_cpm | Label Projection |     6799.78 |           377.203  |    140.556  | zebrafish_labels |
|  0 | accuracy | 0.271506 | logistic_regression_scran   | Label Projection |    14368.6  |            38.7461 |    991.807  | zebrafish_labels |
|  1 | f1       | 0.236048 | mlp_scran                   | Label Projection |    13842.6  |           988.41   |   1126.67   | zebrafish_labels |
|  1 | f1       | 0.247458 | mlp_log_cpm                 | Label Projection |     6325.73 |           376.273  |     89.4715 | zebrafish_labels |
|  1 | f1       | 0.29952  | logistic_regression_log_cpm | Label Projection |     6799.78 |           377.203  |    140.556  | zebrafish_labels |
|  1 | f1       | 0.325879 | logistic_regression_scran   | Label Projection |    14368.6  |            38.7461 |    991.807  | zebrafish_labels |
|  0 | accuracy | 0.824734 | mlp_log_cpm                 | Label Projection |     6204.43 |           376.273  |     52.361  | zebrafish_random |
|  0 | accuracy | 0.829432 | logistic_regression_log_cpm | Label Projection |     6625.18 |           376.273  |     58.9772 | zebrafish_random |
|  0 | accuracy | 0.836516 | mlp_scran                   | Label Projection |    13849    |           610.406  |    898.129  | zebrafish_random |
|  0 | accuracy | 0.841984 | logistic_regression_scran   | Label Projection |    14269.9  |           143.969  |    556.141  | zebrafish_random |
|  1 | f1       | 0.824972 | mlp_log_cpm                 | Label Projection |     6204.43 |           376.273  |     52.361  | zebrafish_random |
|  1 | f1       | 0.826931 | logistic_regression_log_cpm | Label Projection |     6625.18 |           376.273  |     58.9772 | zebrafish_random |
|  1 | f1       | 0.835683 | mlp_scran                   | Label Projection |    13849    |           610.406  |    898.129  | zebrafish_random |
|  1 | f1       | 0.840092 | logistic_regression_scran   | Label Projection |    14269.9  |           143.969  |    556.141  | zebrafish_random |