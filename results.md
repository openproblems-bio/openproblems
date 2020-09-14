|    | metric   |       value | method                               | task                        |   memory_mb |   memory_leaked_mb |   runtime_s | dataset             |
|---:|:---------|------------:|:-------------------------------------|:----------------------------|------------:|-------------------:|------------:|:--------------------|
|  0 | accuracy | 0.96673     | mlp_scran                            | Label Projection            |    11480.3  |         1583.9     |   8032.8    | pancreas_batch      |
|  0 | accuracy | 0.962402    | logistic_regression_scran            | Label Projection            |     9715.58 |         2859.99    |   6131.92   | pancreas_batch      |
|  0 | accuracy | 0.96132     | mlp_log_cpm                          | Label Projection            |    10374.1  |         1193.23    |   1247.76   | pancreas_batch      |
|  0 | accuracy | 0.957533    | logistic_regression_log_cpm          | Label Projection            |     7388.76 |         1234.49    |   3832.58   | pancreas_batch      |
|  1 | f1       | 0.967599    | mlp_scran                            | Label Projection            |    11480.3  |         1583.9     |   8032.8    | pancreas_batch      |
|  1 | f1       | 0.963136    | logistic_regression_scran            | Label Projection            |     9715.58 |         2859.99    |   6131.92   | pancreas_batch      |
|  1 | f1       | 0.961816    | mlp_log_cpm                          | Label Projection            |    10374.1  |         1193.23    |   1247.76   | pancreas_batch      |
|  1 | f1       | 0.957188    | logistic_regression_log_cpm          | Label Projection            |     7388.76 |         1234.49    |   3832.58   | pancreas_batch      |
|  0 | accuracy | 0.988654    | mlp_scran                            | Label Projection            |    11596    |         -129.852   |   2042.24   | pancreas_random     |
|  0 | accuracy | 0.98804     | mlp_log_cpm                          | Label Projection            |    11859.4  |         1192.83    |   1125.12   | pancreas_random     |
|  0 | accuracy | 0.985281    | logistic_regression_scran            | Label Projection            |    11829.9  |         -366.125   |   2839.29   | pancreas_random     |
|  0 | accuracy | 0.984667    | logistic_regression_log_cpm          | Label Projection            |    12093.3  |         1193.82    |   1216.17   | pancreas_random     |
|  1 | f1       | 0.988574    | mlp_scran                            | Label Projection            |    11596    |         -129.852   |   2042.24   | pancreas_random     |
|  1 | f1       | 0.987946    | mlp_log_cpm                          | Label Projection            |    11859.4  |         1192.83    |   1125.12   | pancreas_random     |
|  1 | f1       | 0.985162    | logistic_regression_scran            | Label Projection            |    11829.9  |         -366.125   |   2839.29   | pancreas_random     |
|  1 | f1       | 0.984502    | logistic_regression_log_cpm          | Label Projection            |    12093.3  |         1193.82    |   1216.17   | pancreas_random     |
|  0 | accuracy | 0.250487    | logistic_regression_scran            | Label Projection            |    14136.6  |          267.969   |   2515.71   | zebrafish_labels    |
|  0 | accuracy | 0.241144    | logistic_regression_log_cpm          | Label Projection            |     6784.05 |          377.191   |    476.718  | zebrafish_labels    |
|  0 | accuracy | 0.186259    | mlp_log_cpm                          | Label Projection            |     6371.76 |          376.219   |     96.3792 | zebrafish_labels    |
|  0 | accuracy | 0.182886    | mlp_scran                            | Label Projection            |    13841.3  |          898.965   |   1319.6    | zebrafish_labels    |
|  1 | f1       | 0.296262    | logistic_regression_scran            | Label Projection            |    14136.6  |          267.969   |   2515.71   | zebrafish_labels    |
|  1 | f1       | 0.274748    | logistic_regression_log_cpm          | Label Projection            |     6784.05 |          377.191   |    476.718  | zebrafish_labels    |
|  1 | f1       | 0.199748    | mlp_log_cpm                          | Label Projection            |     6371.76 |          376.219   |     96.3792 | zebrafish_labels    |
|  1 | f1       | 0.197855    | mlp_scran                            | Label Projection            |    13841.3  |          898.965   |   1319.6    | zebrafish_labels    |
|  0 | accuracy | 0.841637    | logistic_regression_scran            | Label Projection            |    14185.7  |          301.129   |   1237.9    | zebrafish_random    |
|  0 | accuracy | 0.838929    | mlp_scran                            | Label Projection            |    13922.1  |          603.562   |   1101.95   | zebrafish_random    |
|  0 | accuracy | 0.82833     | logistic_regression_log_cpm          | Label Projection            |     6540.18 |          376.273   |     50.3583 | zebrafish_random    |
|  0 | accuracy | 0.819975    | mlp_log_cpm                          | Label Projection            |     6183.93 |          376.273   |     64.0222 | zebrafish_random    |
|  1 | f1       | 0.839594    | logistic_regression_scran            | Label Projection            |    14185.7  |          301.129   |   1237.9    | zebrafish_random    |
|  1 | f1       | 0.838426    | mlp_scran                            | Label Projection            |    13922.1  |          603.562   |   1101.95   | zebrafish_random    |
|  1 | f1       | 0.825164    | logistic_regression_log_cpm          | Label Projection            |     6540.18 |          376.273   |     50.3583 | zebrafish_random    |
|  1 | f1       | 0.819633    | mlp_log_cpm                          | Label Projection            |     6183.93 |          376.273   |     64.0222 | zebrafish_random    |
|  0 | knn_auc  | 0.0192371   | mnn_log_scran_pooling                | Multimodal Data Integration |   392598    |          469.57    |    299.433  | scicar_cell_lines   |
|  0 | knn_auc  | 0.00400918  | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    10403    |          933.652   |    355.739  | scicar_cell_lines   |
|  0 | knn_auc  | 0.00269102  | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     7808.71 |          688.152   |    143.401  | scicar_cell_lines   |
|  0 | knn_auc  | 0.00183253  | procrustes                           | Multimodal Data Integration |     5452.34 |            0.34375 |     17.5809 | scicar_cell_lines   |
|  0 | knn_auc  | 0.00158035  | mnn_log_cpm                          | Multimodal Data Integration |   358310    |          134.586   |    101.913  | scicar_cell_lines   |
|  1 | mse      | 1.06542     | mnn_log_scran_pooling                | Multimodal Data Integration |   392598    |          469.57    |    299.433  | scicar_cell_lines   |
|  1 | mse      | 1.00051     | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     7808.71 |          688.152   |    143.401  | scicar_cell_lines   |
|  1 | mse      | 1.00051     | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    10403    |          933.652   |    355.739  | scicar_cell_lines   |
|  1 | mse      | 0.999424    | mnn_log_cpm                          | Multimodal Data Integration |   358310    |          134.586   |    101.913  | scicar_cell_lines   |
|  1 | mse      | 0.890782    | procrustes                           | Multimodal Data Integration |     5452.34 |            0.34375 |     17.5809 | scicar_cell_lines   |
|  0 | knn_auc  | 0.019729    | mnn_log_scran_pooling                | Multimodal Data Integration |   429002    |          736.094   |    473.523  | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000440322 | procrustes                           | Multimodal Data Integration |     6056.96 |           30.9297  |     37.9718 | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000314393 | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |    22236.3  |         3897.51    |   1416.68   | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000131428 | mnn_log_cpm                          | Multimodal Data Integration |   372573    |          136.637   |    274.892  | scicar_mouse_kidney |
|  0 | knn_auc  | 9.31453e-05 | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    22347.9  |         3370.47    |   1796.5    | scicar_mouse_kidney |
|  1 | mse      | 1.04866     | mnn_log_cpm                          | Multimodal Data Integration |   372573    |          136.637   |    274.892  | scicar_mouse_kidney |
|  1 | mse      | 1.04403     | mnn_log_scran_pooling                | Multimodal Data Integration |   429002    |          736.094   |    473.523  | scicar_mouse_kidney |
|  1 | mse      | 1.00022     | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    22347.9  |         3370.47    |   1796.5    | scicar_mouse_kidney |
|  1 | mse      | 1.00022     | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |    22236.3  |         3897.51    |   1416.68   | scicar_mouse_kidney |
|  1 | mse      | 0.942629    | procrustes                           | Multimodal Data Integration |     6056.96 |           30.9297  |     37.9718 | scicar_mouse_kidney |