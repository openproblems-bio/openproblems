|    | metric   |      value | method     | task                        | dataset             |
|---:|:---------|-----------:|:-----------|:----------------------------|:--------------------|
|  0 | accuracy | 0.505859   | dummy      | Label Projection            | dummy               |
|  0 | accuracy | 1          | cheat      | Label Projection            | dummy               |
|  1 | dummy    | 0          | cheat      | Label Projection            | dummy               |
|  1 | dummy    | 0          | dummy      | Label Projection            | dummy               |
|  0 | knn_auc  | 0.00298303 | procrustes | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.0587961  | cheat      | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.062768   | mnn        | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 0          | cheat      | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 0.833141   | procrustes | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 1.05657    | mnn        | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.00119093 | mnn        | Multimodal Data Integration | scicar_mouse_kidney |
|  0 | knn_auc  | 0.0223434  | procrustes | Multimodal Data Integration | scicar_mouse_kidney |
|  0 | knn_auc  | 0.0424163  | cheat      | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0          | cheat      | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0.939055   | procrustes | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0.980936   | mnn        | Multimodal Data Integration | scicar_mouse_kidney |