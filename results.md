|    | metric   |      value | method                               | task                        |   memory_mb |   memory_leaked_mb |   runtime_s | dataset             |
|---:|:---------|-----------:|:-------------------------------------|:----------------------------|------------:|-------------------:|------------:|:--------------------|
|  0 | accuracy | 0.621951   | logistic_regression_log_cpm          | Label Projection            |     7.26562 |         1.48828    |   0.0520863 | pancreas_batch      |
|  0 | accuracy | 0.579268   | mlp_log_cpm                          | Label Projection            |     6.24609 |         0.285156   |   0.579627  | pancreas_batch      |
|  0 | accuracy | 0.567073   | mlp_scran                            | Label Projection            |    30.9023  |        24.9648     |   0.586668  | pancreas_batch      |
|  0 | accuracy | 0.554878   | logistic_regression_scran            | Label Projection            |   633.137   |       627.602      |  14.9918    | pancreas_batch      |
|  1 | f1       | 0.624757   | logistic_regression_log_cpm          | Label Projection            |     7.26562 |         1.48828    |   0.0520863 | pancreas_batch      |
|  1 | f1       | 0.583483   | mlp_log_cpm                          | Label Projection            |     6.24609 |         0.285156   |   0.579627  | pancreas_batch      |
|  1 | f1       | 0.572338   | mlp_scran                            | Label Projection            |    30.9023  |        24.9648     |   0.586668  | pancreas_batch      |
|  1 | f1       | 0.560702   | logistic_regression_scran            | Label Projection            |   633.137   |       627.602      |  14.9918    | pancreas_batch      |
|  0 | accuracy | 0.8        | mlp_log_cpm                          | Label Projection            |     5.96875 |         0          |   0.620226  | pancreas_random     |
|  0 | accuracy | 0.746667   | logistic_regression_log_cpm          | Label Projection            |     5.98047 |         0.0546875  |   0.0506659 | pancreas_random     |
|  0 | accuracy | 0.706667   | mlp_scran                            | Label Projection            |    26.4219  |        20.4375     |   0.626644  | pancreas_random     |
|  0 | accuracy | 0.693333   | logistic_regression_scran            | Label Projection            |     6.13672 |         0.117188   |   0.0792254 | pancreas_random     |
|  1 | f1       | 0.801547   | mlp_log_cpm                          | Label Projection            |     5.96875 |         0          |   0.620226  | pancreas_random     |
|  1 | f1       | 0.749658   | logistic_regression_log_cpm          | Label Projection            |     5.98047 |         0.0546875  |   0.0506659 | pancreas_random     |
|  1 | f1       | 0.710184   | mlp_scran                            | Label Projection            |    26.4219  |        20.4375     |   0.626644  | pancreas_random     |
|  1 | f1       | 0.697304   | logistic_regression_scran            | Label Projection            |     6.13672 |         0.117188   |   0.0792254 | pancreas_random     |
|  0 | accuracy | 0.0163934  | mlp_log_cpm                          | Label Projection            |     5.76562 |         0          |   0.244471  | zebrafish_labels    |
|  0 | accuracy | 0          | mlp_scran                            | Label Projection            |     5.74219 |         0          |   0.211888  | zebrafish_labels    |
|  0 | accuracy | 0          | logistic_regression_scran            | Label Projection            |     5.79688 |         0.0703125  |   0.288389  | zebrafish_labels    |
|  0 | accuracy | 0          | logistic_regression_log_cpm          | Label Projection            |     5.83984 |         0.0976562  |   0.0701691 | zebrafish_labels    |
|  1 | f1       | 0.00252207 | mlp_log_cpm                          | Label Projection            |     5.76562 |         0          |   0.244471  | zebrafish_labels    |
|  1 | f1       | 0          | logistic_regression_scran            | Label Projection            |     5.79688 |         0.0703125  |   0.288389  | zebrafish_labels    |
|  1 | f1       | 0          | logistic_regression_log_cpm          | Label Projection            |     5.83984 |         0.0976562  |   0.0701691 | zebrafish_labels    |
|  1 | f1       | 0          | mlp_scran                            | Label Projection            |     5.74219 |         0          |   0.211888  | zebrafish_labels    |
|  0 | accuracy | 0.192308   | logistic_regression_scran            | Label Projection            |     5.66406 |         0          |   0.0699402 | zebrafish_random    |
|  0 | accuracy | 0.134615   | logistic_regression_log_cpm          | Label Projection            |     5.71484 |         0          |   0.0692928 | zebrafish_random    |
|  0 | accuracy | 0.115385   | mlp_scran                            | Label Projection            |     5.62109 |         0          |   0.22847   | zebrafish_random    |
|  0 | accuracy | 0.115385   | mlp_log_cpm                          | Label Projection            |     5.63672 |         0          |   0.208209  | zebrafish_random    |
|  1 | f1       | 0.154773   | logistic_regression_scran            | Label Projection            |     5.66406 |         0          |   0.0699402 | zebrafish_random    |
|  1 | f1       | 0.136126   | logistic_regression_log_cpm          | Label Projection            |     5.71484 |         0          |   0.0692928 | zebrafish_random    |
|  1 | f1       | 0.107166   | mlp_log_cpm                          | Label Projection            |     5.63672 |         0          |   0.208209  | zebrafish_random    |
|  1 | f1       | 0.0779093  | mlp_scran                            | Label Projection            |     5.62109 |         0          |   0.22847   | zebrafish_random    |
|  0 | knn_auc  | 0.769706   | procrustes                           | Multimodal Data Integration |     5.87891 |         0.015625   |   0.0487141 | scicar_cell_lines   |
|  0 | knn_auc  | 0.627982   | mnn_log_scran_pooling                | Multimodal Data Integration |    10.2695  |         0.332031   |   3.88503   | scicar_cell_lines   |
|  0 | knn_auc  | 0.52371    | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.72656 |         0.015625   |   0.214223  | scicar_cell_lines   |
|  0 | knn_auc  | 0.521865   | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     7.18359 |         1.53125    |   0.216146  | scicar_cell_lines   |
|  0 | knn_auc  | 0.417506   | mnn_log_cpm                          | Multimodal Data Integration |    35.4883  |         3.46875    |   7.68056   | scicar_cell_lines   |
|  1 | mse      | 1.01319    | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.72656 |         0.015625   |   0.214223  | scicar_cell_lines   |
|  1 | mse      | 0.998209   | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     7.18359 |         1.53125    |   0.216146  | scicar_cell_lines   |
|  1 | mse      | 0.984945   | mnn_log_cpm                          | Multimodal Data Integration |    35.4883  |         3.46875    |   7.68056   | scicar_cell_lines   |
|  1 | mse      | 0.974549   | mnn_log_scran_pooling                | Multimodal Data Integration |    10.2695  |         0.332031   |   3.88503   | scicar_cell_lines   |
|  1 | mse      | 0.442196   | procrustes                           | Multimodal Data Integration |     5.87891 |         0.015625   |   0.0487141 | scicar_cell_lines   |
|  0 | knn_auc  | 0.846586   | procrustes                           | Multimodal Data Integration |     6.05859 |        -0.00390625 |   0.0877079 | scicar_mouse_kidney |
|  0 | knn_auc  | 0.59586    | mnn_log_scran_pooling                | Multimodal Data Integration |    14.1914  |         0.0703125  |   3.90184   | scicar_mouse_kidney |
|  0 | knn_auc  | 0.552441   | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.85156 |         0.00390625 |   0.2584    | scicar_mouse_kidney |
|  0 | knn_auc  | 0.549842   | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     6.13672 |         0.132812   |   0.250555  | scicar_mouse_kidney |
|  0 | knn_auc  | 0.41964    | mnn_log_cpm                          | Multimodal Data Integration |    14.1328  |         0.414062   |   3.90117   | scicar_mouse_kidney |
|  1 | mse      | 1.03358    | mnn_log_scran_pooling                | Multimodal Data Integration |    14.1914  |         0.0703125  |   3.90184   | scicar_mouse_kidney |
|  1 | mse      | 1.00711    | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     6.13672 |         0.132812   |   0.250555  | scicar_mouse_kidney |
|  1 | mse      | 1.004      | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.85156 |         0.00390625 |   0.2584    | scicar_mouse_kidney |
|  1 | mse      | 1.00367    | mnn_log_cpm                          | Multimodal Data Integration |    14.1328  |         0.414062   |   3.90117   | scicar_mouse_kidney |
|  1 | mse      | 0.185888   | procrustes                           | Multimodal Data Integration |     6.05859 |        -0.00390625 |   0.0877079 | scicar_mouse_kidney |