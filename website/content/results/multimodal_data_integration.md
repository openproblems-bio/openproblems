+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = "Multimodal Data Integration"
subtitle = "Realigning multimodal measurements of the same cell using various technologies"

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
|    | metric   |       value | method                               | task                        |   memory_mb |   memory_leaked_mb |   runtime_s | dataset             |
|---:|:---------|------------:|:-------------------------------------|:----------------------------|------------:|-------------------:|------------:|:--------------------|
|  0 | knn_auc  | 0.00158035  | mnn_log_cpm                          | Multimodal Data Integration |   358310    |          134.586   |    101.913  | scicar_cell_lines   |
|  0 | knn_auc  | 0.00183253  | procrustes                           | Multimodal Data Integration |     5452.34 |            0.34375 |     17.5809 | scicar_cell_lines   |
|  0 | knn_auc  | 0.00269102  | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     7808.71 |          688.152   |    143.401  | scicar_cell_lines   |
|  0 | knn_auc  | 0.00400918  | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    10403    |          933.652   |    355.739  | scicar_cell_lines   |
|  0 | knn_auc  | 0.0192371   | mnn_log_scran_pooling                | Multimodal Data Integration |   392598    |          469.57    |    299.433  | scicar_cell_lines   |
|  1 | mse      | 0.890782    | procrustes                           | Multimodal Data Integration |     5452.34 |            0.34375 |     17.5809 | scicar_cell_lines   |
|  1 | mse      | 0.999424    | mnn_log_cpm                          | Multimodal Data Integration |   358310    |          134.586   |    101.913  | scicar_cell_lines   |
|  1 | mse      | 1.00051     | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    10403    |          933.652   |    355.739  | scicar_cell_lines   |
|  1 | mse      | 1.00051     | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     7808.71 |          688.152   |    143.401  | scicar_cell_lines   |
|  1 | mse      | 1.06542     | mnn_log_scran_pooling                | Multimodal Data Integration |   392598    |          469.57    |    299.433  | scicar_cell_lines   |
|  0 | knn_auc  | 9.31453e-05 | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    22347.9  |         3370.47    |   1796.5    | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000131428 | mnn_log_cpm                          | Multimodal Data Integration |   372573    |          136.637   |    274.892  | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000314393 | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |    22236.3  |         3897.51    |   1416.68   | scicar_mouse_kidney |
|  0 | knn_auc  | 0.000440322 | procrustes                           | Multimodal Data Integration |     6056.96 |           30.9297  |     37.9718 | scicar_mouse_kidney |
|  0 | knn_auc  | 0.019729    | mnn_log_scran_pooling                | Multimodal Data Integration |   429002    |          736.094   |    473.523  | scicar_mouse_kidney |
|  1 | mse      | 0.942629    | procrustes                           | Multimodal Data Integration |     6056.96 |           30.9297  |     37.9718 | scicar_mouse_kidney |
|  1 | mse      | 1.00022     | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |    22236.3  |         3897.51    |   1416.68   | scicar_mouse_kidney |
|  1 | mse      | 1.00022     | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |    22347.9  |         3370.47    |   1796.5    | scicar_mouse_kidney |
|  1 | mse      | 1.04403     | mnn_log_scran_pooling                | Multimodal Data Integration |   429002    |          736.094   |    473.523  | scicar_mouse_kidney |
|  1 | mse      | 1.04866     | mnn_log_cpm                          | Multimodal Data Integration |   372573    |          136.637   |    274.892  | scicar_mouse_kidney |