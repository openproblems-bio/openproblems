|    | metric   |      value | method     | task                        |
|---:|:---------|-----------:|:-----------|:----------------------------|
|  0 | accuracy | 1          | cheat      | Label Projection            |
|  1 | dummy    | 0          | cheat      | Label Projection            |
|  0 | accuracy | 0.590909   | dummy      | Label Projection            |
|  1 | dummy    | 0          | dummy      | Label Projection            |
|  0 | knn_auc  | 0.697669   | cheat      | Multimodal Data Integration |
|  1 | mse      | 0          | cheat      | Multimodal Data Integration |
|  0 | knn_auc  | 0.00745504 | mnn        | Multimodal Data Integration |
|  1 | mse      | 1.08376    | mnn        | Multimodal Data Integration |
|  0 | knn_auc  | 0.00670292 | procrustes | Multimodal Data Integration |
|  1 | mse      | 0.999936   | procrustes | Multimodal Data Integration |